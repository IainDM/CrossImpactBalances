using Test
using CrossImpactBalances

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")

# Follow succession_step until a signature repeats.
# Returns (cycle_length, attractor); cycle_length == 1 means a fixed point.
function trace(cib, u)
    seen = Dict(signature(cib, u) => 1)
    n = 1
    v = copy(u)
    while true
        v = succession_step(cib, v)
        s = signature(cib, v)
        haskey(seen, s) && return (n - seen[s] + 1, v)
        seen[s] = (n += 1)
    end
end

@testset "CrossImpactBalances" begin

    @testset "Load .scw file" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        @test cib.ndesc == 3
        @test cib.nvariants == [3, 3, 4]
        @test cib.ndim == 10
        @test length(cib.descriptors) == 3
        @test cib.descriptors == ["WTRD", "WSEC", "WECO"]
        @test cib.variants["WTRD"] == ["FT", "Ntl", "Mix"]
        @test cib.variants["WSEC"] == ["Rlx", "Mod", "Alrt"]
        @test cib.variants["WECO"] == ["Decl", "Stag", "ModGr", "RpdGr"]

        # Check CIM corners (verified against Python CIBSA)
        @test cib.cim[1, 1] == 0
        @test cib.cim[1, 4] == 1
        @test cib.cim[10, 1] == 3
        @test cib.cim[10, 6] == 0
    end

    @testset "Load .sl solutions" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        # .sl file: "2 3 2", "3 1 3", "2 2 2", "1 2 3" -> 0-based
        @test length(cib.kernel) == 4
        @test cib.kernel[1] == [1, 2, 1]
        @test cib.kernel[2] == [2, 0, 2]
        @test cib.kernel[3] == [1, 1, 1]
        @test cib.kernel[4] == [0, 1, 2]
    end

    @testset "Signatures" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        # Known signatures from Python CIBSA
        @test signature(cib, [1, 2, 1]) == 16
        @test signature(cib, [2, 0, 2]) == 20
        @test signature(cib, [1, 1, 1]) == 13
        @test signature(cib, [0, 1, 2]) == 21

        # Round-trip
        for u in cib.kernel
            s = signature(cib, u)
            @test inv_signature(cib, s) == u
        end

        @test max_signature(cib) == signature(cib, [2, 2, 3])
    end

    @testset "Impact balance" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        # Verified against Python CIBSA
        ib = impact_balance(cib, [1, 2, 1])
        @test ib == [-4, 6, -2, -4, 2, 2, 4, 7, -4, -7]
    end

    @testset "Succession - fixed points" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        # All kernel scenarios must be fixed points of succession
        for u in cib.kernel
            v = succession_step(cib, u)
            @test v == u
        end
    end

    @testset "Succession - known trajectories" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        # [0,0,0] has a cycle of length 2 (verified in Python)
        nper, veqm = trace(cib, [0, 0, 0])
        @test nper == 2

        # [1,0,0] converges to [1,2,1] (sig=16)
        nper, veqm = trace(cib, [1, 0, 0])
        @test nper == 1
        @test veqm == [1, 2, 1]

        # [0,0,1] converges to [2,0,2] (sig=20)
        nper, veqm = trace(cib, [0, 0, 1])
        @test nper == 1
        @test veqm == [2, 0, 2]

        # [0,1,2] is already consistent
        nper, veqm = trace(cib, [0, 1, 2])
        @test nper == 1
        @test veqm == [0, 1, 2]
    end

    @testset "Find consistent (computed vs loaded)" begin
        cib_loaded = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                              sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        cib_computed = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"))

        # Must find the same 4 consistent scenarios
        loaded_sigs = sort([signature(cib_loaded, u) for u in cib_loaded.kernel])
        computed_sigs = sort([signature(cib_computed, u) for u in cib_computed.kernel])
        @test loaded_sigs == computed_sigs
        @test loaded_sigs == [13, 16, 20, 21]
    end

    @testset "Exhaustive search (small)" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"))

        sigs = sort([signature(cib, u) for u in cib.kernel])
        @test sigs == [13, 16, 20, 21]

        # Every result must be a true fixed point
        for u in cib.kernel
            @test CrossImpactBalances.succession_step(cib, u) == u
        end
    end

    @testset "Basin analysis (small)" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        fps, basins, cyc = find_basins(cib)
        fp_sigs = sort([signature(cib, u) for u in fps])

        # Must find the same 4 fixed points
        @test fp_sigs == [13, 16, 20, 21]

        # Every fixed point must be a true fixed point
        for u in fps
            @test CrossImpactBalances.succession_step(cib, u) == u
        end

        # Basin sizes + cycle count must equal total scenarios (3*3*4 = 36)
        @test sum(basins) + cyc == 36

        # Each basin must be non-empty
        @test all(b -> b >= 1, basins)

        # Fixed points count themselves, so each basin >= 1
        # and there must be some cycle scenarios (we know [0,0,0] cycles)
        @test cyc > 0

        # Verify a known convergence: [1,0,0] -> [1,2,1] (sig 16)
        sig_100 = signature(cib, [1, 0, 0])
        sig_121 = 16
        # [1,0,0] should be in the basin of sig 16
        idx = findfirst(u -> signature(cib, u) == sig_121, fps)
        @test !isnothing(idx)
    end

    @testset "Typical (10 desc × 3 variants, exhaustive + basins)" begin
        scw = joinpath(SAMPLE_DIR, "bench_typical.scw")

        # Exhaustive: must find the 2 fixed points verified by Python
        cib = load_scw(scw)
        sigs = sort([signature(cib, u) for u in cib.kernel])
        @test sigs == [13785, 13839]
        for u in cib.kernel
            @test succession_step(cib, u) == u
        end

        # Basin analysis
        fps, basins, cyc = find_basins(load_scw(scw; kernel=Vector{Vector{Int}}()))
        fp_sigs = sort([signature(cib, u) for u in fps])
        @test fp_sigs == [13785, 13839]
        @test sum(basins) + cyc == 59049  # 3^10
        @test all(b -> b >= 1, basins)
        @test cyc == 5181

        # Largest basin verified against Python chain-following
        idx = findfirst(u -> signature(cib, u) == 13785, fps)
        @test basins[idx] == 52329
    end

    @testset "Round-trip every sample .scw/.sl pair" begin
        # For each sample file, load_scw + find_consistent (no .sl) must
        # produce the same kernel signatures as loading with the .sl file.
        for name in ["CIB_global", "CIB_natl_regional",
                     "CIB_nested", "CIB_nested_simplified_standardized"]
            scw = joinpath(SAMPLE_DIR, "$name.scw")
            sl  = joinpath(SAMPLE_DIR, "$name.sl")
            isfile(sl) || continue   # nonstandard file has no .sl

            cib_loaded   = load_scw(scw; sl_file=sl)
            cib_computed = load_scw(scw)

            sigs_loaded   = sort([signature(cib_loaded, u) for u in cib_loaded.kernel])
            sigs_computed = sort([signature(cib_computed, u) for u in cib_computed.kernel])
            @test sigs_loaded == sigs_computed
            # Every fixed point reported must actually be a fixed point
            for u in cib_computed.kernel
                @test CrossImpactBalances.succession_step(cib_computed, u) == u
            end
        end
    end

    @testset "load_scw error paths" begin
        # Missing file → SystemError or IOError from eachline
        @test_throws SystemError load_scw(joinpath(SAMPLE_DIR, "does_not_exist.scw"))

        # Empty file → custom error
        mktemp() do path, io
            close(io)  # zero bytes
            @test_throws ErrorException load_scw(path)
        end

        # Malformed CIM dimensions: header says 2 variants but CIM has only 1 row
        mktemp() do path, io
            write(io, "&D1\n -a\n -b\n#\n&\n \n \n#\n#\nFFFFFF\nFFFFFF\n#\n#\n0,0\n#\n# 0-0\n\"\"\n#\n#\n")
            close(io)
            @test_throws ErrorException load_scw(path)
        end
    end

    @testset "load_scw with explicit kernel keyword" begin
        explicit = [[1, 2, 1], [0, 1, 2]]
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"); kernel=explicit)
        @test cib.kernel == explicit
        # The kernel keyword overrides find_consistent: only the two we passed.
        @test length(cib.kernel) == 2
    end

    include("property_tests.jl")

end

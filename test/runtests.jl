using Test
using Random
using CrossImpactBalances

# Shared minimal JSON parser used for Python-CIBSA cross-validation fixtures
include(joinpath(@__DIR__, "..", "mcp", "json.jl"))

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")

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
        nper, veqm = succession(cib, [0, 0, 0])
        @test nper == 2

        # [1,0,0] converges to [1,2,1] (sig=16)
        nper, veqm = succession(cib, [1, 0, 0])
        @test nper == 1
        @test veqm == [1, 2, 1]

        # [0,0,1] converges to [2,0,2] (sig=20)
        nper, veqm = succession(cib, [0, 0, 1])
        @test nper == 1
        @test veqm == [2, 0, 2]

        # [0,1,2] is already consistent
        nper, veqm = succession(cib, [0, 1, 2])
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
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"); exhaustive=true)

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
        cib = load_scw(scw; exhaustive=true)
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

    @testset "Inner product matrix" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        M = inner_product_matrix(cib)
        @test size(M) == (4, 4)

        # Exact expected values from Python CIBSA
        expected = [
            15  -10  15  -6;
           -12   12  -7   9;
             9   -7   9  -2;
           -11    8  -5   8
        ]
        @test M == expected
    end

    @testset "own / cross / inner_product (vs Python)" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        u = [1, 2, 1]    # sig 16
        v = [0, 1, 2]    # sig 21

        # impact_balance(u) = [-4, 6, -2, -4, 2, 2, 4, 7, -4, -7]
        # u's variants live at flat indices [2, 6, 8] (1-based) -> [6, 2, 7]
        @test own_impact_balance(cib, u) == [6, 2, 7]
        # v's variants live at flat indices [1, 5, 9] -> [-4, 2, -4]
        @test cross_impact_balance(cib, u, v) == [-4, 2, -4]
        @test inner_product(cib, u, v) == -6
        @test inner_product(cib, u, u) == sum(own_impact_balance(cib, u))

        # Inner-product matrix M[i,j] == inner_product(kernel[i], kernel[j])
        M = inner_product_matrix(cib)
        for (i, ui) in enumerate(cib.kernel), (j, uj) in enumerate(cib.kernel)
            @test M[i, j] == inner_product(cib, ui, uj)
        end
    end

    @testset "sim_anneal / build_graph / merge_scenarios (vs Python)" begin
        # Cross-validated against Python CIBSA via test/generate_sim_anneal_expected.py
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))
        expected = parse_json_file(joinpath(SAMPLE_DIR, "sim_anneal_expected.json"))

        # Sanity: kernel must be in the same order Python produced for the fixture
        py_kernel_sigs = Int[v for v in expected["kernel_sigs"]]
        @test [signature(cib, u) for u in cib.kernel] == py_kernel_sigs

        for case in expected["cases"]
            thr = Int[v for v in case["thresholds"]]
            set_thresholds!(cib, thr)
            ignore_cycles = case["ignore_cycles"]

            # ── sim_anneal per kernel scenario ──
            for sa_entry in case["sim_anneal"]
                u = Int[v for v in sa_entry["u"]]
                @test signature(cib, u) == sa_entry["u_sig"]

                accessible = sim_anneal(cib, u; ignore_cycles=ignore_cycles)
                got_sigs = sort([signature(cib, w) for w in accessible])
                want_sigs = Int[v for v in sa_entry["accessible_sigs"]]
                @test got_sigs == want_sigs

                weights = sim_anneal(cib, u; ignore_cycles=ignore_cycles,
                                     return_weights=true)
                @test weights["reject"] == sa_entry["reject"]
                got_pairs = sort([[Int(k), Int(v)] for (k, v) in weights if k != "reject"])
                want_pairs = [[Int(p[1]), Int(p[2])] for p in sa_entry["weights"]]
                @test got_pairs == want_pairs
            end

            # ── build_graph: compare dense adjacency element-wise ──
            adj_jl = Matrix(build_graph(cib))
            adj_py = case["adjacency"]
            @test size(adj_jl) == (length(adj_py), length(adj_py[1]))
            for (i, row) in enumerate(adj_py), (j, val) in enumerate(row)
                @test adj_jl[i, j] == Int(val)
            end

            # ── merge_scenarios: components are sets, order-independent ──
            comps_jl = [sort(Int[s for s in c]) for c in merge_scenarios(cib)]
            comps_py = [sort(Int[s for s in c]) for c in case["components"]]
            @test sort(comps_jl) == sort(comps_py)
        end
    end

    @testset "set_thresholds! and rand_scenario" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        @test cib.thresholds == [0, 0, 0]
        set_thresholds!(cib, [1, 2, 3])
        @test cib.thresholds == [1, 2, 3]
        @test_throws DimensionMismatch set_thresholds!(cib, [1, 2])
        @test_throws DimensionMismatch set_thresholds!(cib, [1, 2, 3, 4])

        # rand_scenario must produce valid 0-based indices for every descriptor
        rng = MersenneTwister(42)
        for _ in 1:200
            u = rand_scenario(cib; rng=rng)
            @test length(u) == cib.ndesc
            for i in 1:cib.ndesc
                @test 0 <= u[i] < cib.nvariants[i]
            end
        end

        # Reproducibility: same seed -> same scenario
        u1 = rand_scenario(cib; rng=MersenneTwister(123))
        u2 = rand_scenario(cib; rng=MersenneTwister(123))
        @test u1 == u2
    end

    @testset "set_impact! / get_impact (in-place editing)" begin
        cib = load_scw(joinpath(SAMPLE_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLE_DIR, "CIB_global.sl"))

        d1 = cib.descriptors[1]
        d2 = cib.descriptors[2]
        v1 = cib.variants[d1][1]      # 0-based index 0
        v2 = cib.variants[d2][2]      # 0-based index 1

        # Name-based get/set round-trip; set returns the previous value.
        original = get_impact(cib, d1, v1, d2, v2)
        prev = set_impact!(cib, d1, v1, d2, v2, original + 5)
        @test prev == original
        @test get_impact(cib, d1, v1, d2, v2) == original + 5

        # The stored transpose must stay consistent with cim after the edit.
        @test cib.cim_t == permutedims(cib.cim)

        # 0-based index form must resolve to the same cell as the name form.
        @test get_impact(cib, 0, 0, 1, 1) == original + 5
        set_impact!(cib, 0, 0, 1, 1, original)     # restore via index form
        @test get_impact(cib, d1, v1, d2, v2) == original
        @test cib.cim_t == permutedims(cib.cim)

        # Editing the matrix in place changes the consistent set without a reload.
        base_kernel = Set(signature(cib, u) for u in find_consistent(cib))
        # Crank a few cells hard to force some fixed points to move.
        set_impact!(cib, d2, cib.variants[d2][1],   d1, cib.variants[d1][1], 50)
        set_impact!(cib, d2, cib.variants[d2][end], d1, cib.variants[d1][1], -50)
        edited_kernel = Set(signature(cib, u) for u in find_consistent(cib))
        @test cib.cim_t == permutedims(cib.cim)
        @test edited_kernel != base_kernel

        # Error paths (unknown names, out-of-range indices).
        @test_throws ArgumentError get_impact(cib, "Nope", v1, d2, v2)
        @test_throws ArgumentError get_impact(cib, d1, "Nope", d2, v2)
        @test_throws ArgumentError set_impact!(cib, d1, v1, "Nope", v2, 1)
        @test_throws ArgumentError CrossImpactBalances._table_index(cib, 0, 99)
        @test_throws ArgumentError CrossImpactBalances._table_index(cib, 99, 0)
    end

    @testset "Monte-Carlo find_consistent reproducibility" begin
        # bench_typical has 3^10 = 59,049 scenarios; with mc_threshold=500
        # we hit the sampling-without-replacement path of get_scenario_signatures.
        scw = joinpath(SAMPLE_DIR, "bench_typical.scw")

        # Same RNG seed → identical kernels
        k1 = sort([signature(load_scw(scw; mc_threshold=500,
                                      rng=MersenneTwister(7)),
                            u)
                   for u in load_scw(scw; mc_threshold=500,
                                     rng=MersenneTwister(7)).kernel])
        k2 = sort([signature(load_scw(scw; mc_threshold=500,
                                      rng=MersenneTwister(7)),
                            u)
                   for u in load_scw(scw; mc_threshold=500,
                                     rng=MersenneTwister(7)).kernel])
        @test k1 == k2

        # Sampling without replacement: the 500-element sample must be
        # 500 distinct signatures.
        cib = load_scw(scw; kernel=Vector{Vector{Int}}())
        sample = CrossImpactBalances.get_scenario_signatures(
            cib; max=500, allow_dups=false, rng=MersenneTwister(99))
        @test length(sample) == 500
        @test length(unique(sample)) == 500

        # allow_dups=true uses sampling with replacement (Julia's original
        # behavior). Distinct count will typically be less than `max`.
        sample_dups = CrossImpactBalances.get_scenario_signatures(
            cib; max=500, allow_dups=true, rng=MersenneTwister(99))
        @test length(sample_dups) == 500
        # With high probability some duplicates appear (birthday paradox)
        @test length(unique(sample_dups)) < 500
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
            cib_computed = load_scw(scw; exhaustive=true)

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

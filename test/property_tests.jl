# Property tests: cross-check the optimized search/basin implementations
# against brute-force oracles built only on the public scoring primitives
# (succession_step / succession). Instances are constructed directly, so
# these tests cover shapes the .scw fixtures don't: radix-1 descriptors,
# ties-only matrices, all-zero matrices, and entries large enough to force
# the wide-integer arithmetic path.

using Test
using Random
using CrossImpactBalances

# Construct a CIB directly without going through a .scw file.
function make_cib(nvariants::Vector{Int}, cim::Matrix{Int})
    ndesc = length(nvariants)
    ndim = sum(nvariants)
    @assert size(cim) == (ndim, ndim)
    descriptors = ["D$i" for i in 1:ndesc]
    variants = Dict(descriptors[i] => ["V$(i)_$j" for j in 1:nvariants[i]]
                    for i in 1:ndesc)
    desc_offsets = cumsum(vcat(0, nvariants[1:end-1]))
    return CIB(descriptors, variants, nvariants, cim, permutedims(cim),
               ndim, ndesc, Vector{Vector{Int}}(), zeros(Int, ndesc),
               10^9, desc_offsets)
end

# Oracle 1: a scenario is consistent iff succession_step maps it to itself.
function naive_kernel(cib::CIB)
    kern = Vector{Vector{Int}}()
    for s in 0:max_signature(cib)
        u = inv_signature(cib, s)
        succession_step(cib, u) == u && push!(kern, u)
    end
    return kern
end

# Oracle 2: follow succession from every scenario; tally fixed-point basins
# and count starts that end in non-fixed-point cycles.
function naive_basins(cib::CIB)
    tally = Dict{Int,Int}()
    cyc = 0
    for s in 0:max_signature(cib)
        nper, veqm = succession(cib, inv_signature(cib, s))
        if nper > 1
            cyc += 1
        else
            k = signature(cib, veqm)
            tally[k] = get(tally, k, 0) + 1
        end
    end
    fp_sigs = sort!(collect(keys(tally)))
    return fp_sigs, [tally[k] for k in fp_sigs], cyc
end

# Random CIM with the standard zero within-descriptor blocks (a descriptor
# does not impact itself); set zero_diag=false for adversarial instances.
function rand_cim(rng::AbstractRNG, nvariants::Vector{Int};
                  lo::Int=-3, hi::Int=3, zero_diag::Bool=true)
    ndim = sum(nvariants)
    cim = rand(rng, lo:hi, ndim, ndim)
    if zero_diag
        off = 0
        for nv in nvariants
            cim[off+1:off+nv, off+1:off+nv] .= 0
            off += nv
        end
    end
    return cim
end

sorted_sigs(cib, kern) = sort!([signature(cib, u) for u in kern])

# One instance: every implementation must agree with the oracles.
function check_case(cib::CIB)
    want = sorted_sigs(cib, naive_kernel(cib))

    # Threaded exhaustive sweep: exact same set, ascending-signature order.
    got_exh = find_consistent(cib; exhaustive=true)
    got_exh_sigs = [signature(cib, u) for u in got_exh]
    @test got_exh_sigs == want          # sorted == pins the ordering guarantee

    # Succession-walk full enumeration (mc_threshold above space size).
    got_walk = find_consistent(cib; exhaustive=false)
    @test sorted_sigs(cib, got_walk) == want

    # Basin analysis vs chain-following oracle.
    want_fps, want_sizes, want_cyc = naive_basins(cib)
    fps, sizes, cyc = find_basins(cib)
    got = sort!(collect(zip([signature(cib, u) for u in fps], sizes)))
    @test first.(got) == want_fps
    @test last.(got) == want_sizes
    @test cyc == want_cyc
    @test sum(sizes) + cyc == max_signature(cib) + 1
end

@testset "Property: implementations vs naive oracles" begin
    @testset "seeded random instances" begin
        rng = MersenneTwister(20260721)
        for case in 1:50
            ndesc = rand(rng, 2:6)
            nvariants = [rand(rng, 1:5) for _ in 1:ndesc]  # radix-1 allowed
            while prod(nvariants) > 4000                    # keep oracles fast
                nvariants[argmax(nvariants)] -= 1
            end
            cim = rand_cim(rng, nvariants; zero_diag=(case % 5 != 0))
            check_case(make_cib(nvariants, cim))
        end
    end

    @testset "all-zero matrix (everything is a fixed point)" begin
        nvariants = [2, 3, 2]
        cib = make_cib(nvariants, zeros(Int, 7, 7))
        @test length(naive_kernel(cib)) == 12
        check_case(cib)
    end

    @testset "ties-only entries (-1:1)" begin
        rng = MersenneTwister(99)
        for _ in 1:5
            nvariants = [rand(rng, 2:4) for _ in 1:rand(rng, 3:5)]
            cim = rand_cim(rng, nvariants; lo=-1, hi=1)
            check_case(make_cib(nvariants, cim))
        end
    end

    @testset "large entries (forces wide-integer path)" begin
        rng = MersenneTwister(7)
        for _ in 1:3
            nvariants = [rand(rng, 2:3) for _ in 1:3]
            cim = rand_cim(rng, nvariants; lo=-5000, hi=5000)
            check_case(make_cib(nvariants, cim))
        end
    end

    @testset "radix-1 heavy" begin
        rng = MersenneTwister(1234)
        for nvariants in ([1, 3, 2], [2, 1, 1, 3], [1, 1, 4], [3, 1, 2, 1, 2])
            cim = rand_cim(rng, nvariants)
            check_case(make_cib(nvariants, cim))
        end
    end
end

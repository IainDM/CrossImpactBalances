# Property tests: cross-check the optimized search/basin implementations
# against brute-force oracles built only on the public scoring primitive
# succession_step. Instances are constructed directly, so
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
               ndim, ndesc, Vector{Vector{Int}}(), desc_offsets)
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

# Follow succession_step from `u` until a signature repeats. Returns
# (cycle_length, attractor); cycle_length == 1 means a fixed point.
function naive_trace(cib::CIB, u::Vector{Int})
    seen = Dict{Int,Int}(signature(cib, u) => 1)
    history = 1
    v = copy(u)
    while true
        v = succession_step(cib, v)
        s = signature(cib, v)
        haskey(seen, s) && return (history - seen[s] + 1, v)
        history += 1
        seen[s] = history
    end
end

# Oracle 2: follow succession from every scenario; tally fixed-point basins
# and count starts that end in non-fixed-point cycles.
function naive_basins(cib::CIB)
    tally = Dict{Int,Int}()
    cyc = 0
    for s in 0:max_signature(cib)
        nper, veqm = naive_trace(cib, inv_signature(cib, s))
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
    got_exh = find_consistent(cib; algorithm=:sweep)
    got_exh_sigs = [signature(cib, u) for u in got_exh]
    @test got_exh_sigs == want          # sorted == pins the ordering guarantee

    # Branch-and-bound: identical kernel, same ascending order.
    got_bnb = find_consistent(cib; algorithm=:bnb)
    @test [signature(cib, u) for u in got_bnb] == want

    # Tiny node budget: B&B trips (when the tree exceeds one charge batch)
    # and the sweep fallback must deliver the identical kernel either way.
    got_fb = find_consistent(cib; algorithm=:bnb, bnb_node_budget=1)
    @test [signature(cib, u) for u in got_fb] == want

    # Independent oracle: trace from every start and collect the fixed points
    # it lands on. Built only on succession_step, so it shares no code with
    # the sweep or the branch-and-bound search.
    got_walk = Set{Int}()
    for s in 0:max_signature(cib)
        nper, veqm = naive_trace(cib, inv_signature(cib, s))
        nper == 1 && push!(got_walk, signature(cib, veqm))
    end
    @test sort!(collect(got_walk)) == want

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

# ─── Pluggable succession rules ──────────────────────────────────────────────
#
# A test-only rule whose step is exactly global succession, but as a distinct
# type it is forced down the generic (non-specialised) find_consistent /
# find_basins path. Comparing it to the fast GlobalSuccession path proves the
# generic machinery yields identical results to the hand-optimised one.
struct RefGlobal <: SuccessionRule end
CrossImpactBalances.succession_step(::RefGlobal, cib::CIB, u::Vector{Int}) =
    succession_step(GlobalSuccession(), cib, u)

# Independent brute-force oracle for SEQUENTIAL succession — written from
# scratch (Gauss–Seidel: update descriptors one at a time using the impact
# balance of the partially-updated scenario). Validates both halves of the
# rule's contract: its margin-0 claim (kernel == the fast searches' output)
# and its genuinely different dynamics (basins via the generic path).
function seq_step_oracle(cib::CIB, u::Vector{Int})
    v = copy(u)
    for i in 1:cib.numberOfDescriptors
        ib = impact_balance(cib, v)
        off = cib.desc_offsets[i]
        nv = cib.numberOfVariants[i]
        best = v[i]
        bestval = ib[off + v[i] + 1]
        for j in 0:nv-1
            if ib[off + j + 1] > bestval
                bestval = ib[off + j + 1]
                best = j
            end
        end
        v[i] = best
    end
    return v
end

function seq_basins_oracle(cib::CIB)
    tally = Dict{Int,Int}()
    cyc = 0
    for s0 in 0:max_signature(cib)
        seen = Set{Int}()
        v = inv_signature(cib, s0)
        sig = s0
        while true
            push!(seen, sig)
            v2 = seq_step_oracle(cib, v)
            sig2 = signature(cib, v2)
            if sig2 == sig                 # fixed point
                tally[sig2] = get(tally, sig2, 0) + 1
                break
            elseif sig2 in seen            # entered a cycle of length > 1
                cyc += 1
                break
            end
            v = v2
            sig = sig2
        end
    end
    fp = sort!(collect(keys(tally)))
    return fp, [tally[k] for k in fp], cyc
end

# A test-only *threshold* (inertial) rule: a descriptor moves only when some
# variant beats the current one by more than `m`. It declares a
# fixed_point_margin, so find_consistent must route it through the fast
# threaded sweep / branch-and-bound — validated below against a from-scratch
# brute-force scan. `m == 0` must reproduce global succession exactly.
struct Threshold <: SuccessionRule
    m::Int
end
function CrossImpactBalances.succession_step(rule::Threshold, cib::CIB, u::Vector{Int})
    ib = impact_balance(cib, u)
    v = copy(u)
    for i in 1:cib.numberOfDescriptors
        off = cib.desc_offsets[i]
        cur = ib[off + u[i] + 1]
        best = u[i]
        bestval = cur
        for j in 0:cib.numberOfVariants[i]-1
            if ib[off + j + 1] > bestval && ib[off + j + 1] - cur > rule.m
                bestval = ib[off + j + 1]
                best = j
            end
        end
        v[i] = best
    end
    return v
end
CrossImpactBalances.fixed_point_margin(rule::Threshold) = rule.m

# Brute-force fixed points of the threshold rule: every scenario equal to its
# own step. Shares no code with the sweep / branch-and-bound under test.
threshold_fixed_points(cib, m) =
    sort!([sig for sig in 0:max_signature(cib)
           if succession_step(Threshold(m), cib, inv_signature(cib, sig)) ==
              inv_signature(cib, sig)])

@testset "Property: pluggable succession rules" begin
    rng = MersenneTwister(20260722)
    cases = CIB[]
    for _ in 1:20
        ndesc = rand(rng, 2:5)
        nvariants = [rand(rng, 1:4) for _ in 1:ndesc]
        while prod(nvariants) > 2000
            nvariants[argmax(nvariants)] -= 1
        end
        push!(cases, make_cib(nvariants, rand_cim(rng, nvariants)))
    end

    @testset "generic path reproduces the fast GlobalSuccession path" begin
        for cib in cases
            fast = [signature(cib, u) for u in find_consistent(cib)]
            gen  = [signature(cib, u) for u in
                    find_consistent(cib; rule=RefGlobal())]
            @test gen == fast

            f1, s1, c1 = find_basins(cib)
            f2, s2, c2 = find_basins(cib; rule=RefGlobal())
            @test [signature(cib, u) for u in f2] == [signature(cib, u) for u in f1]
            @test s2 == s1
            @test c2 == c1
        end
    end

    @testset "SequentialSuccession vs independent oracle" begin
        differs = false
        for cib in cases
            # Kernel: sequential fixed points coincide with global ones (the
            # induction argument in the SequentialSuccession docstring), so the
            # rule declares fixed_point_margin = 0 and every fast strategy must
            # find exactly the oracle's fixed points, in ascending order.
            want_k = sort!([signature(cib, inv_signature(cib, s))
                            for s in 0:max_signature(cib)
                            if seq_step_oracle(cib, inv_signature(cib, s)) ==
                               inv_signature(cib, s)])
            for alg in (:auto, :sweep, :bnb)
                got_k = [signature(cib, u) for u in
                         find_consistent(cib; rule=SequentialSuccession(),
                                         algorithm=alg)]
                @test got_k == want_k
            end

            # Basins: generic walk matches the from-scratch sequential oracle.
            of, os, oc = seq_basins_oracle(cib)
            gf, gs, gc = find_basins(cib; rule=SequentialSuccession())
            @test [signature(cib, u) for u in gf] == of
            @test gs == os
            @test gc == oc
            @test sum(gs) + gc == max_signature(cib) + 1

            # Confirm the rule is genuinely distinct from global succession on
            # at least one instance (guards against an accidental no-op rule).
            for s in 0:max_signature(cib)
                u = inv_signature(cib, s)
                if succession_step(SequentialSuccession(), cib, u) !=
                   succession_step(GlobalSuccession(), cib, u)
                    differs = true
                    break
                end
            end
        end
        @test differs
    end

    @testset "search-strategy kwargs rejected for rules without a margin" begin
        # RefGlobal declares no fixed_point_margin (it returns nothing), so it
        # has no fast sweep / branch-and-bound and :bnb / :sweep must be
        # refused. (SequentialSuccession no longer serves here: its fixed
        # points provably coincide with global's, so it declares margin 0 and
        # accepts the fast strategies.)
        cib = cases[1]
        @test_throws ArgumentError find_consistent(cib; rule=RefGlobal(),
                                                   algorithm=:bnb)
        @test_throws ArgumentError find_consistent(cib; rule=RefGlobal(),
                                                   algorithm=:sweep)
    end
end

@testset "Property: threshold rules use the fast sweep / branch-and-bound" begin
    rng = MersenneTwister(20260724)
    saw_strict_superset = false
    for _ in 1:25
        ndesc = rand(rng, 2:5)
        nvariants = [rand(rng, 1:4) for _ in 1:ndesc]
        while prod(nvariants) > 2000
            nvariants[argmax(nvariants)] -= 1
        end
        cib = make_cib(nvariants, rand_cim(rng, nvariants))
        nash = threshold_fixed_points(cib, 0)

        # Threshold(0) is exactly global consistency, via the fast path.
        @test nash == [signature(cib, u) for u in find_consistent(cib)]

        for m in (0, 1, 2, 5)
            want = threshold_fixed_points(cib, m)
            rule = Threshold(m)
            # All three fast strategies find exactly the brute-force fixed
            # points, in ascending-signature order.
            for alg in (:auto, :sweep, :bnb)
                got = [signature(cib, u) for u in
                       find_consistent(cib; rule=rule, algorithm=alg)]
                @test got == want
            end
            # A tiny node budget trips branch-and-bound → sweep fallback, and
            # the margin-parameterised sweep must deliver the identical kernel.
            got_fb = [signature(cib, u) for u in
                      find_consistent(cib; rule=rule, algorithm=:bnb, bnb_node_budget=1)]
            @test got_fb == want
            # A larger margin only ever makes more scenarios consistent.
            @test issubset(Set(nash), Set(want))
            length(want) > length(nash) && (saw_strict_superset = true)
        end
    end
    # The margin genuinely makes extra (non-Nash) scenarios sticky somewhere.
    @test saw_strict_superset
end

@testset "Property: B&B pair-difference bound dominates per-column bounds" begin
    # The prune uses, per (chosen, rival) pair, the suffix-summed extreme of
    # the score DIFFERENCE over each undecided descriptor's shared variant
    # choice. That must never be looser than the decoupled bound it replaced
    # (per-column suffix min of the rival minus suffix max of the chosen),
    # because min_a(x_a - y_a) >= min_a x_a - max_a y_a.
    rng = MersenneTwister(20260730)
    for _ in 1:20
        ndesc = rand(rng, 2:6)
        nvariants = [rand(rng, 1:5) for _ in 1:ndesc]
        cib = make_cib(nvariants,
                       rand_cim(rng, nvariants; zero_diag=rand(rng, Bool)))
        sufDiff, pairOffsets = CrossImpactBalances._bnb_bounds(cib)

        # The pre-tightening per-column suffix extremes, rebuilt naively.
        ndim = cib.numberOfDimensions
        sufmin = zeros(Int, ndesc + 1, ndim)
        sufmax = zeros(Int, ndesc + 1, ndim)
        for l in ndesc:-1:1, c in 1:ndim
            rows = cib.desc_offsets[l]+1 : cib.desc_offsets[l]+nvariants[l]
            sufmin[l, c] = sufmin[l+1, c] + minimum(cib.cim[r, c] for r in rows)
            sufmax[l, c] = sufmax[l+1, c] + maximum(cib.cim[r, c] for r in rows)
        end

        for i in 1:ndesc, chosen in 0:nvariants[i]-1, rival in 0:nvariants[i]-1
            p = pairOffsets[i] + chosen * nvariants[i] + rival + 1
            chosenCol = cib.desc_offsets[i] + chosen + 1
            rivalCol = cib.desc_offsets[i] + rival + 1
            for k in 1:ndesc+1
                @test sufDiff[p, k] >= sufmin[k, rivalCol] - sufmax[k, chosenCol]
                chosen == rival && @test sufDiff[p, k] == 0
            end
        end
    end
end

@testset "B&B node-count regression guard (bench_typical)" begin
    # Guards pruning EFFECTIVENESS, which the kernel-equality tests above are
    # blind to. With the pair-difference bound the search expands 3,933 nodes
    # of bench_typical's 59,049 (the per-column bound needed 7,875). The
    # threshold leaves headroom only for the small thread-count-dependent
    # variation in how root prefixes are charged.
    path = joinpath(@__DIR__, "sample_files", "bench_typical.scw")
    cib = load_scw(path; kernel=Vector{Vector{Int}}())
    sufDiff, pairOffsets = CrossImpactBalances._bnb_bounds(cib)
    kernel, nodes = CrossImpactBalances._bnb_fixed_points(cib, sufDiff, pairOffsets;
                                                          node_budget=typemax(Int))
    @test kernel !== nothing
    @test nodes <= 4_500
end

# Large-scenario-space machinery: exact Int128 counting and its overflow
# guard, the find_basins method dispatch and memory guard, the streaming
# range API, the sampling estimator's statistics and its determinism
# contract, and the package's own random number generator.
#
# The RNG pins below are a CONTRACT, not a convenience: estimate_basins
# promises bit-identical results for a seed on any machine, thread count and
# Julia version, which is only true while the generator's stream never
# changes. If an edit to rng.jl changes these values, that edit breaks every
# published seeded estimate — don't update the pins, revert the change.

using Test
using Random
using CrossImpactBalances
const CIBmod = CrossImpactBalances

const SAMPLES_DIR = joinpath(@__DIR__, "sample_files")

@testset "scenario_count and the max_signature guard" begin
    for name in ["CIB_global", "bench_typical"]
        cib = load_scw(joinpath(SAMPLES_DIR, "$name.scw"); kernel=Vector{Vector{Int}}())
        @test scenario_count(cib) == Int128(max_signature(cib)) + 1
        @test scenario_count(cib) == prod(Int128.(cib.numberOfVariants))
    end

    # A model shaped like the real 6.7e12-scenario matrix (11 quaternary + 13
    # ternary descriptors = 2^22 · 3^13). Signatures fit Int64 comfortably.
    nv_n1 = vcat(fill(4, 11), fill(3, 13))
    big = make_cib(nv_n1, zeros(Int, sum(nv_n1), sum(nv_n1)))
    @test scenario_count(big) == Int128(6_687_075_336_192)
    @test max_signature(big) == 6_687_075_336_191

    # 41 ternary descriptors: 3^41 ≈ 3.6e19 > typemax(Int64). The count stays
    # exact; every signature-DRIVEN entry point refuses, naming alternatives.
    huge = make_cib(fill(3, 41), zeros(Int, 123, 123))
    @test scenario_count(huge) == Int128(3)^41
    @test_throws ArgumentError max_signature(huge)
    @test_throws ArgumentError find_basins(huge)
    # NOT plain find_consistent(huge). Branch-and-bound numbers no scenario, so
    # it now accepts a model this size — and on an all-zero matrix every
    # scenario is consistent and nothing prunes, so calling it here would run
    # for centuries. Only the odometer sweep needs signatures, and that is
    # what is asserted. The >Int64 search itself is covered further down, on a
    # model whose kernel is small and provably known.
    @test_throws ArgumentError find_consistent(huge; algorithm=:sweep)
end

@testset "find_consistent past typemax(Int64)" begin
    # EXACTLY 2^63 scenarios. max_signature ANSWERS here — the largest
    # signature is typemax(Int), so the guard does not fire — and it was the
    # `+ 1` that wrapped, sending a negative count past the "< 100_000" test
    # and into a sweep over an empty range: an empty kernel, no error. Counting
    # in Int128 makes this the ordinary case it always looked like.
    boundary = forcing_chain(63)
    @test scenario_count(boundary) == Int128(2)^63
    @test max_signature(boundary) == typemax(Int)
    @test find_consistent(boundary) == [zeros(Int, 63), ones(Int, 63)]

    # 1.18e21 scenarios: past every Int64 signature, and branch-and-bound does
    # not care. Completeness here is the construction's, checked exhaustively
    # against the oracles at n = 7, 12 and 20 in property_tests.jl; soundness
    # is rechecked below against the public primitive.
    chain = forcing_chain(70)
    @test scenario_count(chain) == Int128(2)^70
    @test_throws ArgumentError max_signature(chain)

    kernel = find_consistent(chain)                    # :auto picks bnb
    @test kernel == [zeros(Int, 70), ones(Int, 70)]
    @test find_consistent(chain; algorithm=:bnb) == kernel
    for u in kernel
        @test CIBmod.succession_step(chain, u) == u
    end

    # ORDER. The all-ones scenario's Int64 signature wraps to -1, so the old
    # sort key put it FIRST; comparing variant choices puts it second, where
    # ascending signature order actually puts it.
    @test signature(chain, kernel[2]) == -1            # why that key is unusable
    @test kernel[1] == zeros(Int, 70)
end

@testset "what still refuses past typemax(Int64), and how it reads" begin
    chain = forcing_chain(70)

    # The sweep walks signatures in order, so it cannot go here — and its
    # message must point at the search that can, not at max_signature's
    # blanket refusal, which now names find_consistent as still working.
    swept = try; find_consistent(chain; algorithm=:sweep); nothing; catch e; e; end
    @test swept isa ArgumentError
    for needle in (":sweep", ":bnb", "typemax(Int64)", "fix_descriptor")
        @test occursin(needle, swept.msg)
    end

    # A rule that declares no fixed_point_margin has only the one-at-a-time
    # scan, which is signature-driven too.
    scanned = try; find_consistent(chain; rule=RefGlobal()); nothing; catch e; e; end
    @test scanned isa ArgumentError
    @test occursin("fixed_point_margin", scanned.msg)

    # Branch-and-bound out of budget has nothing to fall back to at this size.
    # An all-zero matrix prunes nothing, so the budget really does trip.
    flat = make_cib(fill(2, 70), zeros(Int, 140, 140))
    spent = try
        find_consistent(flat; algorithm=:bnb, bnb_node_budget=1); nothing
    catch e; e; end
    @test spent isa ArgumentError
    @test occursin("bnb_node_budget", spent.msg)
end

@testset "kernel order: descriptor comparison == signature order" begin
    # The documented ordering guarantee, unchanged for every model that fits
    # Int64 — checked over every ordered pair of scenarios, including radix-1
    # descriptors (whose digit never differs) and ragged radices.
    for nvariants in ([2, 3, 4], [1, 5, 2, 3], [3, 3, 3, 3], [2, 1, 1, 2])
        cib = make_cib(nvariants, zeros(Int, sum(nvariants), sum(nvariants)))
        scenarios = [inv_signature(cib, s) for s in 0:max_signature(cib)]
        @test all(CIBmod._sig_less(a, b) == (signature(cib, a) < signature(cib, b))
                  for a in scenarios, b in scenarios)
    end
end

@testset "find_basins method dispatch and memory guard" begin
    cib = load_scw(joinpath(SAMPLES_DIR, "CIB_global.scw"); kernel=Vector{Vector{Int}}())
    reference = find_basins(cib)
    @test find_basins(cib; method=:table) == reference
    @test find_basins(cib; method=:stream) == reference
    @test_throws ArgumentError find_basins(cib; method=:tables)
    @test_throws ArgumentError find_basins(cib; signature_range=0:5)              # needs :stream
    @test_throws ArgumentError find_basins(cib; method=:stream, signature_range=0:99)  # out of range
    @test_throws ArgumentError find_basins(cib; method=:stream, signature_range=5:4)   # empty

    # The N1-shaped model: :auto must refuse (the tables would need ~73 TiB;
    # no test machine has that) with a message that names every alternative.
    nv_n1 = vcat(fill(4, 11), fill(3, 13))
    rng = Random.MersenneTwister(31)
    big = make_cib(nv_n1, rand_cim(rng, nv_n1))
    guard = try; find_basins(big); nothing; catch e; e; end
    @test guard isa ArgumentError
    for needle in ("estimate_basins", ":stream", "product_basins", "GiB")
        @test occursin(needle, guard.msg)
    end

    # A model whose TABLE BYTES overflow Int64 (3^39 ≈ 4e18 scenarios): even
    # an explicit method=:table cannot exist and must get the same guidance.
    astronomical = make_cib(fill(3, 39), zeros(Int, 117, 117))
    @test_throws ArgumentError find_basins(astronomical; method=:table)

    # Streaming a RANGE of the N1-shaped space works on any machine: the range
    # invariant replaces the whole-space one.
    fps, sizes, cyc = find_basins(big; method=:stream, signature_range=0:9_999)
    @test sum(sizes) + cyc == 10_000
    fps2, sizes2, cyc2 = find_basins(big; method=:stream, signature_range=0:9_999, cache_bytes=0)
    @test (fps, sizes, cyc) == (fps2, sizes2, cyc2)
end

@testset "find_basins at exactly 2^63 scenarios" begin
    # The one size where max_signature ANSWERS but its `+ 1` wraps: the
    # largest signature is exactly typemax(Int), so the guard does not fire,
    # and the count went negative. :table was shielded by accident (its
    # byte-count guard rejects a model this size first); :stream took the
    # negative count into the memo cache's sizing and died on an InexactError
    # — a crash where a refusal belongs. All three entry points must now
    # refuse in the same words, naming what still works at this scale.
    boundary = forcing_chain(63)
    @test scenario_count(boundary) == Int128(2)^63
    @test max_signature(boundary) == typemax(Int)

    for method in (:auto, :table, :stream)
        refusal = try; find_basins(boundary; method=method); nothing; catch e; e; end
        @test refusal isa ArgumentError
        @test occursin("$(Int128(2)^63) scenarios", refusal.msg)
        for needle in ("estimate_basins", "product_basins")
            @test occursin(needle, refusal.msg)
        end
    end

    # The range API validates against that same count, so it refuses too
    # rather than checking a range against a wrapped bound.
    @test_throws ArgumentError find_basins(boundary; method=:stream, signature_range=0:9)

    # One under the line still runs: the count is the only thing that was
    # wrong, and 2^62 scenarios are as walkable (a range at a time) as any.
    justUnder = forcing_chain(62)
    @test scenario_count(justUnder) == Int128(2)^62
    _, sizes, cycles = find_basins(justUnder; method=:stream, signature_range=0:999)
    @test sum(sizes) + cycles == 1_000
end

@testset "package RNG: pinned streams and unbiased bounded draws" begin
    # THE STREAM-STABILITY CONTRACT (see the file header before touching).
    pins = [
        (UInt64(0),          [0x53175d61490b23df, 0x61da6f3dc380d507,
                              0x5c0fdf91ec9a7bfc, 0x02eebf8c3bbe5e1a]),
        (UInt64(42),         [0xd0764d4f4476689f, 0x519e4174576f3791,
                              0xfbe07cfb0c24ed8c, 0xb37d9f600cd835b8]),
        (UInt64(0xDEADBEEF), [0x0c520eb8fea98ede, 0x2b74a6338b80e0e2,
                              0xbe238770c3795322, 0x5f235f98a244ea97]),
    ]
    for (seed, expected) in pins
        rng = CIBmod._Xoshiro256pp(seed)
        @test [CIBmod._rand64!(rng) for _ in 1:4] == expected
    end
    blockRng = CIBmod._block_rng(0xC1BBA512C1BBA512, 3)
    @test CIBmod._rand64!(blockRng) == 0x5d7bdcb55d4bcdd7

    # Same (seed, block) => same stream; different block => different stream.
    a = CIBmod._block_rng(UInt64(7), 1); b = CIBmod._block_rng(UInt64(7), 1)
    c = CIBmod._block_rng(UInt64(7), 2)
    firstPair = (CIBmod._rand64!(a), CIBmod._rand64!(b))
    @test firstPair[1] == firstPair[2]
    @test CIBmod._rand64!(c) != firstPair[1]

    # Bounded draws stay in range with all mass present (coarse uniformity).
    rng = CIBmod._Xoshiro256pp(UInt64(1))
    for bound in (2, 3, 7)
        draws = [Int(CIBmod._rand_below!(rng, UInt64(bound))) for _ in 1:9_000]
        @test extrema(draws) == (0, bound - 1)
        for value in 0:bound-1
            @test count(==(value), draws) > 9_000 ÷ bound ÷ 2
        end
    end
end

@testset "interval statistics" begin
    # Standard normal quantiles, to more precision than any estimate resolves.
    for (p, want) in ((0.975, 1.9599639845), (0.995, 2.5758293035),
                      (0.9995, 3.2905267315), (0.5, 0.0))
        @test abs(CIBmod._inv_normal_cdf(p) - want) < 1e-6
    end
    @test CIBmod._z_for_confidence(0.95) ≈ 1.9599639845 atol=1e-6
    @test_throws ArgumentError CIBmod._z_for_confidence(1.0)

    # Wilson interval: hand-checked values and the zero-hit upper bound.
    lo, hi = CIBmod._wilson(500, 1000, 1.96)
    @test lo ≈ 0.4691 atol=1e-3
    @test hi ≈ 0.5309 atol=1e-3
    lo0, hi0 = CIBmod._wilson(0, 1000, 1.96)
    @test lo0 == 0.0
    @test 0.0 < hi0 < 0.005          # ≈ z²/n: the honest "no larger than" bound
    @test_throws ArgumentError CIBmod._wilson(5, 0, 1.96)
end

@testset "estimate_basins" begin
    @testset "pinned regression: default seed reproduces exactly" begin
        # These counts define the determinism contract end to end (RNG, walker,
        # block partition, merge). They must be identical at ANY thread count —
        # the CI matrix's 1- and 4-thread legs both run this block.
        cib = load_scw(joinpath(SAMPLES_DIR, "CIB_global.scw"); kernel=Vector{Vector{Int}}())
        est = estimate_basins(cib; samples=50_000)
        @test [signature(cib, u) for u in est.fixedPoints] == [13, 16, 20, 21]
        @test est.hits == [1369, 5483, 16720, 2788]
        @test est.cycleHits == 23640
        @test sum(est.hits) + est.cycleHits == est.samples

        typical = load_scw(joinpath(SAMPLES_DIR, "bench_typical.scw"); kernel=Vector{Vector{Int}}())
        estTypical = estimate_basins(typical; samples=100_000)
        @test [signature(typical, u) for u in estTypical.fixedPoints] == [13785, 13839]
        @test estTypical.hits == [88753, 2619]
        @test estTypical.cycleHits == 8628
    end

    @testset "intervals bracket the exact shares" begin
        typical = load_scw(joinpath(SAMPLES_DIR, "bench_typical.scw"); kernel=Vector{Vector{Int}}())
        fps, sizes, cyc = find_basins(typical)
        n = max_signature(typical) + 1
        est = estimate_basins(typical; samples=200_000, confidence=0.999)
        exactShare = Dict(signature(typical, u) => size / n for (u, size) in zip(fps, sizes))
        for (u, lo, hi) in zip(est.fixedPoints, est.ciLow, est.ciHigh)
            share = exactShare[signature(typical, u)]
            @test lo <= share <= hi
        end
        @test est.cycleCiLow <= cyc / n <= est.cycleCiHigh
        @test est.scenarioCount == scenario_count(typical)
        @test est.sizeEstimates ≈ est.shares .* Float64(n)
    end

    @testset "the kernel is pre-registered: zero-hit attractors stay visible" begin
        # A start biased so heavily that one fixed point's basin is far below
        # 1/samples: it must still appear, with hits == 0 and a positive upper
        # bound. bench_50x50's smallest basin is 3 of 60.5M — the documented
        # CIBSA miss — but too slow for the default suite; CIB_global with few
        # samples shows the same mechanism.
        cib = load_scw(joinpath(SAMPLES_DIR, "CIB_global.scw"),
                       sl_file=joinpath(SAMPLES_DIR, "CIB_global.sl"))
        est = estimate_basins(cib; samples=20)   # signature 13's basin is 1/36
        @test length(est.fixedPoints) == 4        # all four, hit or not
        zeroHit = [i for i in eachindex(est.hits) if est.hits[i] == 0]
        for i in zeroHit
            @test est.ciHigh[i] > 0.0
        end
    end

    @testset "estimation works where exact analysis cannot start" begin
        nv_n1 = vcat(fill(4, 11), fill(3, 13))
        rng = Random.MersenneTwister(31)
        big = make_cib(nv_n1, rand_cim(rng, nv_n1))
        # kernel=[] skips find_consistent (this synthetic model's kernel is not
        # the subject here); 10k samples over 6.7e12 scenarios, in seconds.
        est = estimate_basins(big; samples=10_000, kernel=Vector{Vector{Int}}())
        @test sum(est.hits) + est.cycleHits == 10_000
        @test est.scenarioCount == Int128(6_687_075_336_192)
        # Every attractor sampling found really is a fixed point.
        for u in est.fixedPoints
            @test succession_step(big, u) == u
        end
        # And the printed form renders without error.
        rendered = sprint(show, MIME"text/plain"(), est)
        @test occursin("Exact kernel; estimated shares", rendered)
    end

    @testset "argument validation" begin
        cib = load_scw(joinpath(SAMPLES_DIR, "CIB_global.scw"); kernel=Vector{Vector{Int}}())
        @test_throws ArgumentError estimate_basins(cib; samples=0)
        @test_throws ArgumentError estimate_basins(cib; confidence=1.0)
        @test_throws ArgumentError estimate_basins(cib; kernel=[[0, 0]])        # wrong length
        @test_throws ArgumentError estimate_basins(cib; kernel=[[0, 0, 0]])     # not a fixed point
    end
end

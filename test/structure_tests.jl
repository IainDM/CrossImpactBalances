# Structure analysis: the influence map, island decomposition, dial slicing
# and product composition, held to hand-built cases and to the brute-force
# oracles in property_tests.jl (which is included before this file).

using Test
using Random
using CrossImpactBalances

# A block-diagonal CIM over the given per-island descriptor groups: judgments
# exist only inside an island. `island_of[i]` names descriptor i's island.
function block_diag_cim(rng::AbstractRNG, nvariants::Vector{Int}, island_of::Vector{Int})
    ndim = sum(nvariants)
    offsets = cumsum(vcat(0, nvariants[1:end-1]))
    cim = zeros(Int, ndim, ndim)
    for i in eachindex(nvariants), j in eachindex(nvariants)
        (island_of[i] == island_of[j] && i != j) || continue
        rows = offsets[i]+1:offsets[i]+nvariants[i]
        cols = offsets[j]+1:offsets[j]+nvariants[j]
        cim[rows, cols] .= rand(rng, -3:3, length(rows), length(cols))
    end
    return cim
end

@testset "Influence structure" begin
    @testset "edges exist exactly where variant rows differ" begin
        # D1 (2 variants) votes on D2 (2 variants) with IDENTICAL rows -> no
        # real say; D2 votes on D1 with differing rows -> a real say.
        nv = [2, 2]
        cim = [0 0 5 5;
               0 0 5 5;      # D1's two rows agree on D2's columns (and on D1's)
               1 2 0 0;
               3 4 0 0]      # D2's rows differ on D1's columns
        cib = make_cib(nv, cim)
        st = influence_structure(cib)
        @test st.activeEdges[1, 2] == false     # identical rows: inaudible
        @test st.activeEdges[2, 1] == true
        @test length(st.components) == 1        # the D2->D1 arrow connects them
    end

    @testset "self-influence via a nonzero, non-uniform diagonal block" begin
        nv = [2, 2]
        cim = zeros(Int, 4, 4)
        cim[1, 2] = 7                            # D1 variant 1 promotes variant 2 of itself
        cib = make_cib(nv, cim)
        st = influence_structure(cib)
        @test st.activeEdges[1, 1] == true
        # D1 influences itself, so it is NOT a dial even though no OTHER
        # descriptor has a say over it.
        @test 1 ∉ st.tieFrozen
        @test all(f -> f[1] != 1, st.forced)
    end

    @testset "dials vs jump-and-stays" begin
        # D2 receives no votes at all -> constant level scores -> a dial.
        # D3 receives IDENTICAL rows from D1 that favour its variant 2 (0-based
        # 1) -> constant but unlevel -> jump-and-stay with that one maximizer.
        nv = [2, 3, 2]
        cim = zeros(Int, 7, 7)
        cim[1, 6] = 4; cim[1, 7] = 1             # D1 row 1 votes on D3: [4, 1]
        cim[2, 6] = 4; cim[2, 7] = 1             # D1 row 2 identical -> no say, but a bias
        cim[6, 1] = 1; cim[7, 1] = 2             # D3 has a say over D1 (rows differ)
        cib = make_cib(nv, cim)
        st = influence_structure(cib)
        @test st.tieFrozen == [2]
        @test st.forced == [(3, [0])]
        # Consistency claim of jump-and-stay: every consistent scenario holds a
        # maximizer variant of D3.
        for u in find_consistent(cib)
            @test u[3] == 0
        end
        # Equalising ALL the votes turns forced into a dial (level scores):
        cim2 = copy(cim); cim2[1, 7] = 4; cim2[2, 7] = 4
        st2 = influence_structure(make_cib(nv, cim2))
        @test sort(st2.tieFrozen) == [2, 3]
        @test isempty(st2.forced)

        # Multiple maximizers must all be reported: three variants, two tied
        # on top (identical rows again, so still no real say — just a bias).
        nv3 = [2, 3, 3]
        cim3 = zeros(Int, 8, 8)
        cim3[1, 6] = 4; cim3[1, 7] = 4; cim3[1, 8] = 1
        cim3[2, 6] = 4; cim3[2, 7] = 4; cim3[2, 8] = 1
        cim3[6, 1] = 1; cim3[7, 1] = 2               # D3 has a say over D1
        st3 = influence_structure(make_cib(nv3, cim3))
        @test st3.forced == [(3, [0, 1])]
        @test st3.tieFrozen == [2]
    end

    @testset "dial slicing is exact" begin
        rng = MersenneTwister(20260807)
        for trial in 1:8
            ndesc = rand(rng, 3:5)
            nvariants = [rand(rng, 2:3) for _ in 1:ndesc]
            dial = rand(rng, 1:ndesc)
            cim = rand_cim(rng, nvariants)
            offsets = cumsum(vcat(0, nvariants[1:end-1]))
            # Silence every vote ON the dial descriptor: level constant scores.
            cim[:, offsets[dial]+1:offsets[dial]+nvariants[dial]] .= 0
            cib = make_cib(nvariants, cim)
            st = influence_structure(cib)
            @test dial in st.tieFrozen

            # The pinned models' basins, laid side by side, must reproduce the
            # full model's exactly: per-fixed-point sizes and added cycles.
            want_fps, want_sizes, want_cyc = naive_basins(cib)
            got = Dict{Int,Int}()
            got_cyc = 0
            for variant in 0:nvariants[dial]-1
                fps, sizes, cyc = find_basins(fix_descriptor(cib, dial - 1, variant))
                got_cyc += cyc
                for (u, size) in zip(fps, sizes)
                    full = copy(u); full[dial] = variant
                    sig = signature(cib, full)
                    got[sig] = get(got, sig, 0) + size
                end
            end
            @test sort!(collect(keys(got))) == want_fps
            @test [got[k] for k in sort!(collect(keys(got)))] == want_sizes
            @test got_cyc == want_cyc
        end
    end
end

@testset "split_cib and product_basins" begin
    rng = MersenneTwister(20260808)
    saw_multi_fp = false
    for trial in 1:12
        # Two or three islands of 1-3 descriptors each.
        islands = rand(rng, 2:3)
        island_of = Int[]
        for island in 1:islands
            append!(island_of, fill(island, rand(rng, 1:3)))
        end
        nvariants = [rand(rng, 2:3) for _ in island_of]
        while prod(nvariants) > 3000
            nvariants[argmax(nvariants)] -= 1
        end
        cib = make_cib(nvariants, block_diag_cim(rng, nvariants, island_of))
        st = influence_structure(cib)
        # Islands can only merge if a random block came out all-equal-rows —
        # they can never falsely split. The composition below must be exact
        # for WHATEVER decomposition was found.
        @test length(st.components) >= 1
        parts = split_cib(cib; structure=st)
        @test length(parts) == length(st.components)
        for part in parts
            @test part.numberOfDescriptors == cib.numberOfDescriptors
        end

        composed = product_basins(cib; structure=st)
        want_fps, want_sizes, want_cyc = naive_basins(cib)
        @test [signature(cib, u) for u in composed.fixedPoints] == want_fps
        @test composed.basinSizes == Int128.(want_sizes)
        @test composed.cycleCount == Int128(want_cyc)
        @test sum(composed.basinSizes; init=Int128(0)) + composed.cycleCount ==
              composed.scenarioCount == scenario_count(cib)
        saw_multi_fp |= length(composed.fixedPoints) > 1
    end
    # The trials must have genuinely exercised composition, not just empty or
    # singleton kernels agreeing trivially.
    @test saw_multi_fp

    @testset "custom rules are refused" begin
        cib = make_cib([2, 2], zeros(Int, 4, 4))
        @test_throws ArgumentError product_basins(cib; rule=RefGlobal())
    end

    @testset "structure/model mismatch is refused" begin
        cib3 = make_cib([2, 2, 2], zeros(Int, 6, 6))
        cib2 = make_cib([2, 2], zeros(Int, 4, 4))
        @test_throws ArgumentError split_cib(cib3; structure=influence_structure(cib2))
    end

    @testset "SequentialSuccession composes too" begin
        rng2 = MersenneTwister(9)
        nvariants = [2, 2, 2, 3]
        cib = make_cib(nvariants, block_diag_cim(rng2, nvariants, [1, 1, 2, 2]))
        composed = product_basins(cib; rule=SequentialSuccession())
        of, os, oc = seq_basins_oracle(cib)
        @test [signature(cib, u) for u in composed.fixedPoints] == of
        @test composed.basinSizes == Int128.(os)
        @test composed.cycleCount == Int128(oc)
    end
end

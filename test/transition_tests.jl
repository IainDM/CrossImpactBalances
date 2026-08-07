# The lever map: transition_graph, its world node, and the DOT export.
#
# Every edge this analysis reports is a claim of the form "make this change to
# that consistent scenario and the system ends up here". The tests below check
# those claims the only way worth checking them: by re-walking each one with
# the independent `naive_trace` oracle from property_tests.jl, which shares no
# code with the walker in src/walk.jl.
#
# This file is included after property_tests.jl and structure_tests.jl and
# reuses their helpers (make_cib, rand_cim, naive_trace, seq_step_oracle,
# RefGlobal).

using Test
using Random
using CrossImpactBalances

const TRANSITION_SAMPLES = joinpath(@__DIR__, "sample_files")

# The oracle's view of where a scenario ends up: (isFixedPoint, attractor),
# with a cycle named by its lowest-signature member — the same canonical
# choice src/walk.jl makes, derived here from succession_step alone.
function oracle_destination(cib::CIB, start::Vector{Int})
    cycleLength, attractor = naive_trace(cib, start)
    cycleLength == 1 && return (true, attractor)
    best = attractor
    walker = copy(attractor)
    for _ in 1:cycleLength-1
        walker = succession_step(cib, walker)
        signature(cib, walker) < signature(cib, best) && (best = walker)
    end
    return (false, best)
end

# Apply an edge's commitment to its source scenario.
function apply_changes(graph::TransitionGraph, edge::Transition)
    candidate = copy(graph.nodes[edge.from])
    for (descriptorIndex, variantIndex) in edge.changes
        candidate[descriptorIndex] = variantIndex
    end
    return candidate
end

node_index_of(graph, cib, scenario) =
    findfirst(==(signature(cib, scenario)), [Int(s) for s in graph.signatures])

@testset "Transition graph on CIB_global" begin
    cib = load_scw(joinpath(TRANSITION_SAMPLES, "CIB_global.scw"))
    graph = transition_graph(cib)

    @testset "nodes: the kernel, plus whatever the nudges reached" begin
        attractorSignatures = sort([Int(graph.signatures[i]) for i in 1:length(graph.nodes)
                                    if graph.kinds[i] === :attractor])
        @test attractorSignatures == [13, 16, 20, 21]        # the known kernel
        @test graph.worldIndex == 0                          # no from= was given
        # Signatures ascend, and every node's stored scenario matches its signature.
        @test issorted(graph.signatures)
        for index in 1:length(graph.nodes)
            @test Int(graph.signatures[index]) == signature(cib, graph.nodes[index])
        end
        # Cycle nodes must be the lowest-signature member of their own cycle.
        for index in 1:length(graph.nodes)
            graph.kinds[index] === :cycle || continue
            _, canonical = oracle_destination(cib, graph.nodes[index])
            @test canonical == graph.nodes[index]
        end
    end

    @testset "every attractor × every single nudge, against the oracle" begin
        for sourceIndex in 1:length(graph.nodes)
            graph.kinds[sourceIndex] === :attractor || continue
            isnan(graph.robustness[sourceIndex]) && continue     # destination-only
            source = graph.nodes[sourceIndex]
            expectedReturned = 0
            expectedMoved = 0
            expectedDestabilized = 0
            expectedEdges = Tuple{Int,Int,Int}[]      # (descriptor, variant, destination node)
            for descriptorIndex in 1:cib.numberOfDescriptors,
                variantIndex in 0:cib.numberOfVariants[descriptorIndex]-1
                variantIndex == source[descriptorIndex] && continue
                candidate = copy(source)
                candidate[descriptorIndex] = variantIndex
                isFixedPoint, destination = oracle_destination(cib, candidate)
                if destination == source
                    expectedReturned += 1
                else
                    isFixedPoint ? (expectedMoved += 1) : (expectedDestabilized += 1)
                    push!(expectedEdges, (descriptorIndex, variantIndex,
                                          node_index_of(graph, cib, destination)))
                end
            end
            @test graph.returned[sourceIndex] == expectedReturned
            @test graph.moved[sourceIndex] == expectedMoved
            @test graph.destabilized[sourceIndex] == expectedDestabilized
            # 3+3+4 variants over 3 descriptors = 7 single-descriptor changes.
            @test expectedReturned + expectedMoved + expectedDestabilized == 7
            @test graph.robustness[sourceIndex] ≈ expectedReturned / 7

            got = [(edge.changes[1][1], edge.changes[1][2], edge.to)
                   for edge in graph.edges if edge.from == sourceIndex]
            @test got == expectedEdges          # same edges, same order
        end
    end

    @testset "the current state of the world" begin
        # Trajectories pinned in runtests.jl: [1,0,0] → sig 16, [0,0,1] → sig 20.
        for (world, wantSignature) in ([1, 0, 0] => 16, [0, 0, 1] => 20)
            worldGraph = transition_graph(cib; from=world)
            @test worldGraph.worldIndex == length(worldGraph.nodes)   # world sorts last
            @test worldGraph.kinds[worldGraph.worldIndex] === :world
            @test worldGraph.nodes[worldGraph.worldIndex] == world
            baseline = only([edge for edge in worldGraph.edges
                             if edge.from == worldGraph.worldIndex && isempty(edge.changes)])
            @test Int(worldGraph.signatures[baseline.to]) == wantSignature
            @test worldGraph.kinds[baseline.to] === :attractor
            # The world is a starting point, never a destination.
            @test all(edge -> edge.to != worldGraph.worldIndex, worldGraph.edges)
            @test !isnan(worldGraph.robustness[worldGraph.worldIndex])
        end

        # [0,0,0] cycles (runtests.jl:149) — its default future is a cycle node.
        cyclingGraph = transition_graph(cib; from=[0, 0, 0])
        baseline = only([edge for edge in cyclingGraph.edges
                         if edge.from == cyclingGraph.worldIndex && isempty(edge.changes)])
        @test cyclingGraph.kinds[baseline.to] === :cycle
        _, canonical = oracle_destination(cib, [0, 0, 0])
        @test cyclingGraph.nodes[baseline.to] == canonical
    end

    @testset "a world state that is itself consistent" begin
        # [1,2,1] is consistent (signature 16): it must not become a second
        # node, and it gets no baseline edge to itself.
        graphAt = transition_graph(cib; from=[1, 2, 1])
        @test !any(kind -> kind === :world, graphAt.kinds)   # it is already a future
        @test graphAt.kinds[graphAt.worldIndex] === :attractor
        @test Int(graphAt.signatures[graphAt.worldIndex]) == 16
        @test length(graphAt.nodes) == length(graph.nodes)     # same node set as no-world
        @test all(edge -> !isempty(edge.changes), graphAt.edges)
    end

    @testset "determinism" begin
        @test transition_graph(cib; from=[1, 0, 0]) == transition_graph(cib; from=[1, 0, 0])
        @test transition_graph(cib; radius=2) == transition_graph(cib; radius=2)
    end
end

@testset "Transition graph: from= resolution and argument errors" begin
    cib = load_scw(joinpath(TRANSITION_SAMPLES, "CIB_global.scw"))

    # Dict form (names and 0-based indices mixed) must equal the vector form.
    byName = transition_graph(cib; from=Dict("WTRD" => "Ntl", "WSEC" => "Rlx",
                                             "WECO" => "Decl"))
    byIndex = transition_graph(cib; from=Dict(0 => 1, 1 => 0, 2 => 0))
    @test byName == transition_graph(cib; from=[1, 0, 0])
    @test byIndex == byName

    @test_throws ArgumentError transition_graph(cib; from=Dict("WTRD" => "Ntl"))     # incomplete
    @test_throws ArgumentError transition_graph(cib; from=Dict("WTRD" => "Ntl", 0 => "FT",
                                                               "WSEC" => "Rlx",
                                                               "WECO" => "Decl"))   # twice
    @test_throws ArgumentError transition_graph(cib; from=Dict("Nope" => "Ntl", "WSEC" => "Rlx",
                                                               "WECO" => "Decl"))   # bad name
    @test_throws ArgumentError transition_graph(cib; from=Dict("WTRD" => "Nope", "WSEC" => "Rlx",
                                                               "WECO" => "Decl"))   # bad variant
    @test_throws ArgumentError transition_graph(cib; from=[1, 0])                   # short
    @test_throws ArgumentError transition_graph(cib; from=[1, 0, 9])                # out of range
    @test_throws ArgumentError transition_graph(cib; radius=3)
    @test_throws ArgumentError transition_graph(cib; radius=0)
    # Nothing to analyse at all.
    @test_throws ArgumentError transition_graph(cib; kernel=Vector{Vector{Int}}())

    # An empty kernel plus a world state is a legitimate world-only graph.
    worldOnly = transition_graph(cib; kernel=Vector{Vector{Int}}(), from=[0, 0, 0])
    @test worldOnly.worldIndex != 0
    @test count(kind -> kind === :world, worldOnly.kinds) == 1
    @test !isempty(worldOnly.edges)

    # Duplicate kernel entries collapse to one node.
    duplicated = transition_graph(cib; kernel=[[1, 2, 1], [1, 2, 1], [2, 0, 2]])
    @test count(kind -> kind === :attractor, duplicated.kinds) >= 2
    @test length(unique(duplicated.signatures)) == length(duplicated.signatures)
end

@testset "Transition graph: property tests against the oracles" begin
    rng = MersenneTwister(20260812)
    for trial in 1:12
        ndesc = rand(rng, 2:4)
        nvariants = [rand(rng, 1:3) for _ in 1:ndesc]      # radix-1 allowed
        while prod(nvariants) > 500
            nvariants[argmax(nvariants)] -= 1
        end
        cib = make_cib(nvariants, rand_cim(rng, nvariants; zero_diag=(trial % 4 != 0)))
        kernel = find_consistent(cib)
        world = [rand(rng, 0:count-1) for count in nvariants]
        isempty(kernel) && isempty(world) && continue
        graph = transition_graph(cib; kernel=kernel, from=world, radius=2)

        totalSingles = sum(nvariants) - ndesc

        @testset "trial $trial" begin
            # Node bookkeeping.
            @test length(unique(graph.signatures)) == length(graph.signatures)
            nonWorld = [i for i in 1:length(graph.nodes) if i != graph.worldIndex]
            @test issorted(graph.signatures[nonWorld])
            for index in 1:length(graph.nodes)
                @test Int(graph.signatures[index]) == signature(cib, graph.nodes[index])
                if graph.kinds[index] === :cycle
                    _, canonical = oracle_destination(cib, graph.nodes[index])
                    @test canonical == graph.nodes[index]
                end
            end

            # Every edge re-verified through the independent oracle.
            for edge in graph.edges
                isempty(edge.changes) && continue
                isFixedPoint, destination = oracle_destination(cib, apply_changes(graph, edge))
                @test destination == graph.nodes[edge.to]
                @test graph.kinds[edge.to] === (isFixedPoint ? :attractor : :cycle)
                @test edge.steps >= 1
            end

            # Robustness counts cover exactly the single-descriptor changes.
            for index in 1:length(graph.nodes)
                if isnan(graph.robustness[index])
                    @test graph.returned[index] == graph.moved[index] ==
                          graph.destabilized[index] == 0
                    continue
                end
                @test graph.returned[index] + graph.moved[index] +
                      graph.destabilized[index] == totalSingles
                @test graph.robustness[index] ≈
                      (totalSingles == 0 ? 1.0 : graph.returned[index] / totalSingles)
            end

            # The pair filter, re-derived from scratch: a pair is reported iff
            # it reaches somewhere no single change (nor the baseline) reached.
            for sourceIndex in 1:length(graph.nodes)
                isnan(graph.robustness[sourceIndex]) && continue
                source = graph.nodes[sourceIndex]
                baselineSignature = if graph.kinds[sourceIndex] === :world
                    _, destination = oracle_destination(cib, source)
                    signature(cib, destination)
                else
                    signature(cib, source)
                end
                reachedBySingles = Set{Int}((baselineSignature,))
                for d in 1:ndesc, v in 0:nvariants[d]-1
                    v == source[d] && continue
                    candidate = copy(source); candidate[d] = v
                    _, destination = oracle_destination(cib, candidate)
                    push!(reachedBySingles, signature(cib, destination))
                end
                wantedPairs = Tuple{Int,Int,Int,Int,Int}[]
                for d1 in 1:ndesc-1, v1 in 0:nvariants[d1]-1
                    v1 == source[d1] && continue
                    for d2 in d1+1:ndesc, v2 in 0:nvariants[d2]-1
                        v2 == source[d2] && continue
                        candidate = copy(source); candidate[d1] = v1; candidate[d2] = v2
                        _, destination = oracle_destination(cib, candidate)
                        signature(cib, destination) in reachedBySingles && continue
                        push!(wantedPairs, (d1, v1, d2, v2,
                                            node_index_of(graph, cib, destination)))
                    end
                end
                gotPairs = [(edge.changes[1][1], edge.changes[1][2],
                             edge.changes[2][1], edge.changes[2][2], edge.to)
                            for edge in graph.edges
                            if edge.from == sourceIndex && length(edge.changes) == 2]
                @test gotPairs == wantedPairs
            end
        end
    end
end

@testset "Transition graph: succession rules" begin
    rng = MersenneTwister(20260813)
    nvariants = [3, 2, 3]
    cib = make_cib(nvariants, rand_cim(rng, nvariants))
    kernel = find_consistent(cib)

    # The generic (per-rule) path must reproduce the fast GlobalSuccession path
    # exactly — same nodes, same edges, same step counts.
    @test transition_graph(cib; kernel=kernel, rule=RefGlobal(), from=[0, 0, 0]) ==
          transition_graph(cib; kernel=kernel, from=[0, 0, 0])

    # SequentialSuccession has genuinely different dynamics; check its edges
    # against a from-scratch Gauss–Seidel trace.
    sequential = transition_graph(cib; rule=SequentialSuccession(), from=[0, 0, 0])
    for edge in sequential.edges
        candidate = apply_changes(sequential, edge)
        seen = Set{Int}()
        while !(signature(cib, candidate) in seen)
            push!(seen, signature(cib, candidate))
            next = seq_step_oracle(cib, candidate)
            next == candidate && break
            candidate = next
        end
        if sequential.kinds[edge.to] === :attractor
            @test candidate == sequential.nodes[edge.to]
        else
            @test signature(cib, sequential.nodes[edge.to]) in seen
        end
    end
end

@testset "Transition graph: DOT export" begin
    cib = load_scw(joinpath(TRANSITION_SAMPLES, "CIB_global.scw"))
    dot = to_dot(transition_graph(cib; from=[0, 0, 0]))
    @test occursin("digraph transitions {", dot)
    @test occursin("shape=box", dot)
    @test occursin("style=dashed", dot)               # the cycle node
    @test occursin("current world", dot)
    @test occursin("(default)", dot)                  # the world's baseline edge
    @test occursin("WTRD: ", dot)                     # a commitment label, by name
    @test count(==('{'), dot) == count(==('}'), dot)
    @test endswith(strip(dot), "}")

    # Quotes in names must not break the DOT string literals.
    quoted = CIB(["D\"1", "D2"],
                 Dict("D\"1" => ["a\"a", "b"], "D2" => ["c", "d"]),
                 [2, 2], zeros(Int, 4, 4), zeros(Int, 4, 4), 4, 2,
                 Vector{Vector{Int}}(), [0, 2])
    quotedDot = to_dot(transition_graph(quoted; kernel=[[0, 0]], from=[1, 1]))
    @test occursin("\\\"", quotedDot)
    # Every quote in the output is either an escaped one inside a label or a
    # delimiter — so the unescaped ones must pair up. (An unescaped quote
    # leaking out of a name would leave an odd count and a broken file.)
    unescapedQuotes = 0
    previous = ' '
    for character in quotedDot
        character == '"' && previous != '\\' && (unescapedQuotes += 1)
        # A doubled backslash is an escaped backslash, not an escape prefix.
        previous = (previous == '\\' && character == '\\') ? ' ' : character
    end
    @test iseven(unescapedQuotes)
    # And the escaping lands where it should: a commitment label naming a
    # descriptor and variants that both contain quotes.
    @test occursin("label=\"D\\\"1: b→a\\\"a\"", quotedDot)
end

@testset "Transition graph: cost is independent of scenario-space size" begin
    # The N₁ shape: 11 four-variant + 13 three-variant descriptors = 2^22·3^13
    # ≈ 6.7e12 scenarios. Basin enumeration cannot start here; the lever map
    # does not care, because it only ever walks from known attractors.
    #
    # The kernel is supplied EXPLICITLY: this near-silent matrix has an
    # astronomical number of consistent scenarios, so find_consistent must
    # never be invited to enumerate them.
    nvariants = vcat(fill(4, 11), fill(3, 13))
    ndim = sum(nvariants)
    offsets = cumsum(vcat(0, nvariants[1:end-1]))
    cim = zeros(Int, ndim, ndim)
    cim[1, offsets[2] + 1] = 1                     # one lone influence
    cib = make_cib(nvariants, cim)
    @test scenario_count(cib) == Int128(6_687_075_336_192)

    world = [descriptorIndex % 2 for descriptorIndex in 1:24]
    graph = transition_graph(cib; kernel=[zeros(Int, 24)], from=world, radius=2)

    @test graph.worldIndex == length(graph.nodes)
    @test graph.kinds[graph.worldIndex] === :world
    totalSingles = sum(nvariants) - 24            # 11*3 + 13*2 = 59
    @test totalSingles == 59
    for index in 1:length(graph.nodes)
        isnan(graph.robustness[index]) && continue
        @test graph.returned[index] + graph.moved[index] +
              graph.destabilized[index] == totalSingles
    end
    # Spot-check a handful of edges against the oracle (walking all of them
    # would be the expensive part, not the graph).
    for edge in first(graph.edges, 12)
        isempty(edge.changes) && continue
        isFixedPoint, destination = oracle_destination(cib, apply_changes(graph, edge))
        @test destination == graph.nodes[edge.to]
        @test graph.kinds[edge.to] === (isFixedPoint ? :attractor : :cycle)
    end
    @test to_dot(graph) isa AbstractString
end

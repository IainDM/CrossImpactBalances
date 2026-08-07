# ══ Levers: which commitment moves one future into another ═════════════════
#
# Basin sizes answer "how much of the possibility space drains here" — a
# question about volume. Past a certain model size that answer goes flat:
# hundreds of consistent scenarios, none holding more than a percent or two,
# and the ranking stops being informative. The question that stays sharp is
# about CONNECTIVITY rather than volume:
#
#     We are in this future. What single change puts us in that one?
#
# That is what this file computes. Take a consistent scenario, change one
# descriptor to a different variant, and let succession run (walk.jl). Either
# the system returns to where it started — the commitment did not take, and
# that is a measure of how robust that future is — or it settles somewhere
# else, and you have found a lever: a specific, nameable change that moves the
# system from one consistent future to another.
#
# Do that for every descriptor of every consistent scenario and the result is
# a directed graph whose nodes are the model's futures and whose edges are
# labelled with the commitments that cause the transitions. Add the analyst's
# CURRENT STATE OF THE WORLD as one more node — a scenario that need not be
# consistent — and you also get its default future (where things go if nothing
# is done) and which commitments redirect it.
#
# THE COST DOES NOT DEPEND ON THE SIZE OF THE SCENARIO SPACE. Every walk
# starts from a known attractor and follows the dynamics; nothing enumerates.
# A model with 3.5 quintillion scenarios costs exactly as much as one with
# thirty, because the price is (number of attractors + 1) × (number of
# alternative variants) short walks. That is what makes this the analysis of
# choice where the exact basin analysis cannot run at all.
#
# ── Radius 2 ──
# Some futures are stable against every single change and still reachable by
# changing two things at once. Passing `radius=2` also tries every PAIR of
# descriptors, and — the important part — reports a pair only when it reaches
# somewhere no single change reached. Those are the genuine combination
# levers, the ones invisible to one-at-a-time analysis.

"""
    Transition

One edge of a [`TransitionGraph`](@ref): the commitment `changes`, applied to
node `from`, leads the system to node `to`.

- `changes`: the commitment, as `(descriptor, variant)` pairs — descriptor
  1-based (its position in `cib.descriptors`), variant 0-based (as scenarios
  are numbered everywhere in this package). One pair for a single-descriptor
  lever, two for a pair lever, and **empty for the world node's baseline
  edge**, which is where the current state of the world drifts on its own.
- `steps`: succession steps the walk computed on the way. A rough measure of
  how far the destination was, not a path length — the walker counts its cycle
  canonicalisation lap too (see walk.jl).
"""
struct Transition
    from::Int
    to::Int
    changes::Vector{Tuple{Int,Int}}
    steps::Int
end

"""
    TransitionGraph

The lever map of a model, from [`transition_graph`](@ref).

Nodes (parallel vectors `nodes`, `kinds`, `signatures`, sorted by ascending
signature, with the world node — if any — last):

- `kinds[i]` is `:attractor` for a consistent scenario, `:cycle` for a
  repeating loop reached by some commitment (named by its lowest-signature
  member, so every walk entering that loop names it identically), or `:world`
  for the supplied current state of the world.
- `signatures` are `Int128`, so node identity survives models whose scenario
  count passes `Int64` (see [`scenario_count`](@ref)).

Edges: `edges`, a `Vector{Transition}` ordered by source node, then by the
order the commitments were tried (baseline first, then single descriptors,
then pairs).

Robustness (per node, counting **single-descriptor** commitments only, so the
numbers mean the same thing whatever `radius` was used):

- `returned[i]` — commitments after which the system came back to where it
  started. `moved[i]` — commitments that landed on a different consistent
  scenario. `destabilized[i]` — commitments that landed in a cycle.
- `robustness[i]` = `returned / total`, i.e. the fraction of single changes
  this future absorbs. **`NaN` marks a node that was never perturbed** — a
  cycle, or a fixed point discovered as somebody's destination rather than
  supplied in the kernel. Those nodes are destinations only.

`worldIndex` is the world node's position, or `0` when no `from=` was given.
"""
struct TransitionGraph
    nodes::Vector{Vector{Int}}
    kinds::Vector{Symbol}
    signatures::Vector{Int128}
    edges::Vector{Transition}
    returned::Vector{Int}
    moved::Vector{Int}
    destabilized::Vector{Int}
    robustness::Vector{Float64}
    radius::Int
    worldIndex::Int
    descriptors::Vector{String}
    variantNames::Vector{Vector{String}}
end

# Structural equality, so a determinism check is just `a == b`. (Julia's
# default for mutable-free structs would compare field by field already, but
# spelling it out documents that the graph IS its contents — two runs of the
# same analysis must be indistinguishable.)
Base.:(==)(a::Transition, b::Transition) =
    a.from == b.from && a.to == b.to && a.changes == b.changes && a.steps == b.steps
Base.hash(t::Transition, h::UInt) =
    hash(t.steps, hash(t.changes, hash(t.to, hash(t.from, hash(:Transition, h)))))

function Base.:(==)(a::TransitionGraph, b::TransitionGraph)
    for field in fieldnames(TransitionGraph)
        # `robustness` holds NaN for unperturbed nodes, and NaN != NaN, so
        # compare it with `isequal` (which treats NaNs as equal to each other).
        isequal(getfield(a, field), getfield(b, field)) || return false
    end
    return true
end
Base.hash(g::TransitionGraph, h::UInt) =
    hash(g.edges, hash(g.signatures, hash(g.kinds, hash(:TransitionGraph, h))))

"""
    transition_graph(cib; rule=GlobalSuccession(), kernel=nothing,
                     from=nothing, radius=1) -> TransitionGraph

Map the levers between a model's futures: from every consistent scenario, try
changing each descriptor to each of its other variants, follow succession to
wherever that leads, and record which commitments move the system where.

`from` adds the **current state of the world** as an extra node — a scenario
that need not be consistent. Its first edge is the baseline: where the world
drifts with no commitment at all. The rest are the commitments that redirect
it. Give it either a `Vector{Int}` of 0-based variant indices, one per
descriptor, or a `Dict` naming every descriptor:

    transition_graph(cib; from = Dict("Trade" => "Protectionist",
                                      "Economy" => "Low", ...))

Every descriptor must be assigned — a half-specified world state has no
well-defined trajectory, so it is an error rather than a guess.

`radius=2` additionally tries every pair of descriptors, reporting a pair only
when it reaches a destination no single change reached — the combination
levers. Cost per source node is `Σ(mᵢ-1)` single commitments and about
`(Σ(mᵢ-1))² / 2` pairs, for variant counts `mᵢ`; nothing here scales with the
size of the scenario space, so this analysis runs identically on a model of
thirty scenarios and one of 10¹⁸.

`kernel` supplies the consistent scenarios instead of computing them (they are
taken exactly, as in [`estimate_basins`](@ref); an empty kernel plus a `from`
gives a world-only graph). Use [`to_dot`](@ref) to render the result.
"""
function transition_graph(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                          kernel::Union{Nothing,Vector{Vector{Int}}}=nothing,
                          from::Union{Nothing,AbstractVector{<:Integer},AbstractDict}=nothing,
                          radius::Integer=1)
    radius in (1, 2) || throw(ArgumentError(
        "transition_graph: radius must be 1 (single commitments) or 2 " *
        "(also pairs), got $radius"))
    kernelScenarios, _ = _resolve_kernel(cib, rule, kernel; caller="transition_graph")
    worldScenario = from === nothing ? nothing : _resolve_world(cib, from)
    (isempty(kernelScenarios) && worldScenario === nothing) && throw(ArgumentError(
        "transition_graph: this model has no consistent scenarios to map, and no " *
        "current state of the world was given — pass from=... to analyse a starting " *
        "state, or kernel=... if you know the fixed points"))

    # ── Node registry: signature -> provisional index ──
    # Nodes appear both up front (the kernel, the world) and during the walks
    # (destinations nobody supplied), so the per-node tallies grow with it.
    idForSignature = Dict{Int128,Int}()
    nodeScenarios = Vector{Vector{Int}}()
    nodeKinds = Symbol[]
    nodeSignatures = Int128[]
    returned = Int[]
    moved = Int[]
    destabilized = Int[]
    perturbed = Bool[]
    function register!(scenario::Vector{Int}, kind::Symbol)
        signatureValue = _signature128(cib, scenario)
        existing = get(idForSignature, signatureValue, 0)
        existing == 0 || return existing        # a signature determines its scenario
        push!(nodeScenarios, copy(scenario))
        push!(nodeKinds, kind)
        push!(nodeSignatures, signatureValue)
        push!(returned, 0); push!(moved, 0); push!(destabilized, 0)
        push!(perturbed, false)
        idForSignature[signatureValue] = length(nodeScenarios)
        return length(nodeScenarios)
    end

    # Kernel first (deduplicated), then the world — which may turn out to BE
    # one of the consistent scenarios, in which case it is that same node and
    # keeps its :attractor kind; there is one node per scenario, always.
    sourceIndices = Int[register!(scenario, :attractor) for scenario in kernelScenarios]
    unique!(sourceIndices)
    worldIndex = 0
    if worldScenario !== nothing
        worldIndex = register!(worldScenario, :world)
    end

    # ── Walking apparatus (one buffer set; this loop is not hot enough to thread) ──
    # If radius-2 on a very large kernel ever measures slow, the drop-in is the
    # per-worker _WalkSpace + @sync/@spawn idiom in estimate.jl — with the
    # caveat that the merge must reproduce this loop's edge order exactly.
    scoreType = _score_type(cib)
    cimTranspose = rule isa GlobalSuccession ? Matrix{scoreType}(cib.cim_t) :
                                               Matrix{scoreType}(undef, 0, 0)
    space = _WalkSpace(cib, scoreType)
    walk = _walker(rule, cib, space, cimTranspose)

    edges = Transition[]
    candidate = Vector{Int}(undef, cib.numberOfDescriptors)

    # Sources: the consistent scenarios, then the world if it is a node of its
    # own. Destinations discovered along the way are never perturbed — that is
    # what keeps the cost exactly (kernel + 1) sources however large the space.
    sources = copy(sourceIndices)
    if worldIndex != 0 && nodeKinds[worldIndex] === :world
        push!(sources, worldIndex)
    end

    for sourceIndex in sources
        source = nodeScenarios[sourceIndex]
        perturbed[sourceIndex] = true

        # The baseline: where this node goes with no commitment at all. A
        # consistent scenario goes nowhere by definition (it is its own
        # successor — _resolve_kernel has already verified that), so only the
        # world needs walking, and only the world gets a baseline edge.
        baselineSignature = nodeSignatures[sourceIndex]
        if nodeKinds[sourceIndex] === :world
            isFixedPoint, destination, steps = walk(source)
            baselineIndex = register!(destination, isFixedPoint ? :attractor : :cycle)
            baselineSignature = nodeSignatures[baselineIndex]
            # A self-loop says nothing ("the world stays where it is"), and it
            # is what happens when the world state is itself an attractor or
            # the named member of its own cycle.
            baselineIndex == sourceIndex ||
                push!(edges, Transition(sourceIndex, baselineIndex, Tuple{Int,Int}[], steps))
        end

        # ── Single commitments ──
        # `seen` collects this source's radius-1 destinations, so that the pair
        # pass below can tell a genuinely new destination from one already
        # reachable by changing a single thing.
        seen = Set{Int128}((baselineSignature,))
        for descriptorIndex in 1:cib.numberOfDescriptors
            for variantIndex in 0:cib.numberOfVariants[descriptorIndex]-1
                variantIndex == source[descriptorIndex] && continue
                copyto!(candidate, source)
                candidate[descriptorIndex] = variantIndex
                isFixedPoint, destination, steps = walk(candidate)
                destinationIndex = register!(destination, isFixedPoint ? :attractor : :cycle)
                push!(seen, nodeSignatures[destinationIndex])
                if nodeSignatures[destinationIndex] == baselineSignature
                    returned[sourceIndex] += 1        # the commitment did not take
                else
                    if isFixedPoint
                        moved[sourceIndex] += 1
                    else
                        destabilized[sourceIndex] += 1
                    end
                    push!(edges, Transition(sourceIndex, destinationIndex,
                                            [(descriptorIndex, variantIndex)], steps))
                end
            end
        end

        # ── Pairs, when asked for: only genuinely new destinations ──
        radius == 2 || continue
        for firstDescriptor in 1:cib.numberOfDescriptors-1,
            firstVariant in 0:cib.numberOfVariants[firstDescriptor]-1
            firstVariant == source[firstDescriptor] && continue
            for secondDescriptor in firstDescriptor+1:cib.numberOfDescriptors,
                secondVariant in 0:cib.numberOfVariants[secondDescriptor]-1
                secondVariant == source[secondDescriptor] && continue
                copyto!(candidate, source)
                candidate[firstDescriptor] = firstVariant
                candidate[secondDescriptor] = secondVariant
                isFixedPoint, destination, steps = walk(candidate)
                destinationSignature = _signature128(cib, destination)
                # `seen` is frozen at the radius-1 destinations, so several
                # pairs reaching the same new place are all reported — they are
                # different levers to the same future.
                destinationSignature in seen && continue
                destinationIndex = register!(destination, isFixedPoint ? :attractor : :cycle)
                push!(edges, Transition(sourceIndex, destinationIndex,
                                        [(firstDescriptor, firstVariant),
                                         (secondDescriptor, secondVariant)], steps))
            end
        end
    end

    return _assemble_graph(cib, nodeScenarios, nodeKinds, nodeSignatures, edges,
                           returned, moved, destabilized, perturbed, radius, worldIndex)
end

# One calling convention for both walkers: hand it a scenario, get back
# (isFixedPoint, attractor, steps). The fast path's attractor lives in a reused
# buffer, so it is copied before it escapes; the generic path already returns a
# fresh vector.
function _walker(::GlobalSuccession, cib::CIB, space::_WalkSpace{ScoreInt},
                 cimTranspose::Matrix{ScoreInt}) where {ScoreInt<:Signed}
    return function (start::Vector{Int})
        isFixedPoint, steps = _walk_to_attractor!(space, start, cimTranspose,
                                                  cib.numberOfVariants, cib.desc_offsets,
                                                  cib.numberOfDescriptors,
                                                  cib.numberOfDimensions)
        return (isFixedPoint, copy(space.current), steps)
    end
end

function _walker(rule::SuccessionRule, cib::CIB, ::_WalkSpace, ::Matrix)
    return function (start::Vector{Int})
        return _walk_to_attractor(rule, cib, start)
    end
end

# Put the nodes in their final order — ascending signature, world last — and
# renumber the edges to match. Sorting the non-world nodes keeps the output
# independent of the order things happened to be discovered in.
function _assemble_graph(cib::CIB, nodeScenarios, nodeKinds, nodeSignatures, edges,
                         returned, moved, destabilized, perturbed, radius, worldIndex)
    numberOfNodes = length(nodeScenarios)
    isWorldNode = [index == worldIndex && nodeKinds[index] === :world
                   for index in 1:numberOfNodes]
    order = sort(1:numberOfNodes; by = index -> (isWorldNode[index], nodeSignatures[index]))
    newIndexOf = Vector{Int}(undef, numberOfNodes)
    for (newIndex, oldIndex) in enumerate(order)
        newIndexOf[oldIndex] = newIndex
    end

    totalNudges(scenario) = sum(cib.numberOfVariants) - cib.numberOfDescriptors
    robustness = Vector{Float64}(undef, numberOfNodes)
    for index in 1:numberOfNodes
        if !perturbed[index]
            robustness[index] = NaN            # destination only: never tried
        else
            total = totalNudges(nodeScenarios[index])
            # A model whose descriptors all have a single variant admits no
            # commitments at all; nothing can dislodge such a scenario.
            robustness[index] = total == 0 ? 1.0 : returned[index] / total
        end
    end

    remappedEdges = [Transition(newIndexOf[edge.from], newIndexOf[edge.to],
                                edge.changes, edge.steps) for edge in edges]
    # Edges were emitted source by source in signature order already; sorting
    # by the new source index is a stable no-op in the common case and keeps
    # the promise exact when a world node's edges were emitted last.
    sort!(remappedEdges; by = edge -> edge.from, alg = MergeSort)

    variantNames = [cib.variants[cib.descriptors[descriptorIndex]]
                    for descriptorIndex in 1:cib.numberOfDescriptors]
    return TransitionGraph(nodeScenarios[order], nodeKinds[order], nodeSignatures[order],
                           remappedEdges, returned[order], moved[order],
                           destabilized[order], robustness[order], radius,
                           worldIndex == 0 ? 0 : newIndexOf[worldIndex],
                           copy(cib.descriptors), variantNames)
end

# ── Resolving the current state of the world ────────────────────────────────

"""
    _resolve_world(cib, from) -> Vector{Int}

Turn a user-supplied current state of the world into a scenario (0-based
variant indices). Accepts a plain vector of 0-based indices, or a dictionary
naming descriptors and variants — either by name or by 0-based index, exactly
as [`get_impact`](@ref) does. Every descriptor must appear: unlike
[`fix_descriptor`](@ref), which pins one descriptor and leaves the rest free,
a trajectory has to start from a complete scenario, so a partial world state
is reported as an error rather than filled in with a guess.
"""
function _resolve_world(cib::CIB, from::AbstractVector{<:Integer})
    length(from) == cib.numberOfDescriptors || throw(ArgumentError(
        "transition_graph: from has $(length(from)) entries but the model has " *
        "$(cib.numberOfDescriptors) descriptors — give one variant per descriptor"))
    world = Vector{Int}(undef, cib.numberOfDescriptors)
    for descriptorIndex in 1:cib.numberOfDescriptors
        # _table_index validates the range and names the descriptor if it fails.
        _table_index(cib, descriptorIndex - 1, from[descriptorIndex])
        world[descriptorIndex] = Int(from[descriptorIndex])
    end
    return world
end

function _resolve_world(cib::CIB, from::AbstractDict)
    world = Vector{Int}(undef, cib.numberOfDescriptors)
    assigned = falses(cib.numberOfDescriptors)
    for (descriptorKey, variantKey) in from
        descriptorPosition = _desc_index(cib, descriptorKey)
        assigned[descriptorPosition] && throw(ArgumentError(
            "transition_graph: descriptor \"$(cib.descriptors[descriptorPosition])\" is " *
            "assigned twice in from (by name and by index, or by two equal keys)"))
        # _table_index returns the flat matrix row; subtracting the descriptor's
        # offset turns it back into a 0-based variant number.
        world[descriptorPosition] =
            _table_index(cib, descriptorKey, variantKey) - cib.desc_offsets[descriptorPosition] - 1
        assigned[descriptorPosition] = true
    end
    if !all(assigned)
        missingNames = cib.descriptors[.!assigned]
        throw(ArgumentError(
            "transition_graph: from must assign every descriptor — missing " *
            "$(join(missingNames, ", ")). A current state of the world with gaps has no " *
            "trajectory to follow."))
    end
    return world
end

# ── Rendering ───────────────────────────────────────────────────────────────

# A scenario as readable variant names: "FT/Mod/ModGr".
_scenario_names(graph::TransitionGraph, nodeIndex::Int) =
    join([graph.variantNames[descriptorIndex][graph.nodes[nodeIndex][descriptorIndex] + 1]
          for descriptorIndex in 1:length(graph.descriptors)], "/")

# One commitment as "Trade: FT→Prot", or several joined with " + ".
function _changes_text(graph::TransitionGraph, edge::Transition; separator::AbstractString=" + ")
    isempty(edge.changes) && return "(no commitment — baseline)"
    parts = map(edge.changes) do (descriptorIndex, variantIndex)
        current = graph.nodes[edge.from][descriptorIndex]
        string(graph.descriptors[descriptorIndex], ": ",
               graph.variantNames[descriptorIndex][current + 1], "→",
               graph.variantNames[descriptorIndex][variantIndex + 1])
    end
    return join(parts, separator)
end

function Base.show(io::IO, ::MIME"text/plain", graph::TransitionGraph)
    attractorCount = count(==(:attractor), graph.kinds)
    cycleCount = count(==(:cycle), graph.kinds)
    println(io, "TransitionGraph: ", attractorCount, " consistent scenario(s), ",
            cycleCount, " cycle destination(s)",
            graph.worldIndex == 0 ? "" : ", plus the current state of the world",
            " · radius ", graph.radius, " · ", length(graph.edges), " lever(s)")

    edgesFrom(nodeIndex) = [edge for edge in graph.edges if edge.from == nodeIndex]
    analysed = [index for index in 1:length(graph.nodes)
                if !isnan(graph.robustness[index]) && index != graph.worldIndex]
    # Most robust first: the futures hardest to dislodge are the headline.
    for nodeIndex in sort(analysed; by = index -> -graph.robustness[index])
        total = graph.returned[nodeIndex] + graph.moved[nodeIndex] + graph.destabilized[nodeIndex]
        println(io, "  [", nodeIndex, "] ", _scenario_names(graph, nodeIndex),
                " — robustness ", _pct(graph.robustness[nodeIndex]),
                " (", graph.returned[nodeIndex], "/", total, " single changes absorbed; ",
                graph.moved[nodeIndex], " move it, ", graph.destabilized[nodeIndex],
                " destabilise it)")
        outgoing = edgesFrom(nodeIndex)
        if isempty(outgoing)
            println(io, "        no lever moves this future")
            continue
        end
        for destination in unique(edge.to for edge in outgoing)
            marker = graph.kinds[destination] === :cycle ? " (cycle)" : ""
            println(io, "     → [", destination, "] ", _scenario_names(graph, destination),
                    marker, " via:")
            for edge in outgoing
                edge.to == destination || continue
                println(io, "          ", _changes_text(graph, edge),
                        length(edge.changes) > 1 ? "   [pair]" : "")
            end
        end
    end

    destinationOnly = [index for index in 1:length(graph.nodes) if isnan(graph.robustness[index])]
    if !isempty(destinationOnly)
        println(io, "  Destination-only nodes (not perturbed): ",
                join(["[$index] " * _scenario_names(graph, index) *
                      (graph.kinds[index] === :cycle ? " (cycle)" : "")
                      for index in destinationOnly], ", "))
    end

    if graph.worldIndex != 0
        println(io, "  ── Current state of the world ──")
        println(io, "  [", graph.worldIndex, "] ", _scenario_names(graph, graph.worldIndex))
        if graph.kinds[graph.worldIndex] !== :world
            println(io, "        (this state is itself consistent scenario [",
                    graph.worldIndex, "] — it is already a future, not a starting point)")
        else
            worldEdges = edgesFrom(graph.worldIndex)
            baseline = findfirst(edge -> isempty(edge.changes), worldEdges)
            if baseline === nothing
                println(io, "        default future: stays where it is")
            else
                edge = worldEdges[baseline]
                println(io, "        default future, with no commitment: [", edge.to, "] ",
                        _scenario_names(graph, edge.to),
                        graph.kinds[edge.to] === :cycle ? " (cycle — it never settles)" : "",
                        " after ", edge.steps, " step(s)")
            end
            total = graph.returned[graph.worldIndex] + graph.moved[graph.worldIndex] +
                    graph.destabilized[graph.worldIndex]
            println(io, "        ", graph.returned[graph.worldIndex], "/", total,
                    " single commitments leave that default unchanged; the rest redirect it:")
            for edge in worldEdges
                isempty(edge.changes) && continue
                println(io, "          ", _changes_text(graph, edge), "  ⟶  [", edge.to, "] ",
                        _scenario_names(graph, edge.to),
                        graph.kinds[edge.to] === :cycle ? " (cycle)" : "",
                        length(edge.changes) > 1 ? "   [pair]" : "")
            end
        end
    end
    print(io, "  (to_dot(graph) renders this for Graphviz)")
end

Base.show(io::IO, graph::TransitionGraph) =
    print(io, "TransitionGraph(", length(graph.nodes), " nodes, ", length(graph.edges),
          " levers, radius ", graph.radius, ")")

# ── Graphviz export ─────────────────────────────────────────────────────────

# Quotes and backslashes are the only characters that can break a DOT string
# literal; descriptor and variant names come from user files, so escape them.
_dot_escape(text::AbstractString) =
    replace(replace(String(text), "\\" => "\\\\"), "\"" => "\\\"")

"""
    to_dot(graph) -> String

Render a [`TransitionGraph`](@ref) as a Graphviz `digraph`, ready to paste
into any Graphviz viewer (`dot -Tsvg levers.dot -o levers.svg`, or one of the
browser-based renderers).

Consistent scenarios are boxes labelled with their variant names and
robustness; cycles are dashed ellipses; the current state of the world is a
heavy box. Each edge carries the commitment that causes the transition, so the
picture reads directly as "change this, end up there". Written by hand rather
than through a graph library — the package has no dependencies.
"""
function to_dot(graph::TransitionGraph)
    io = IOBuffer()
    println(io, "digraph transitions {")
    println(io, "  rankdir=LR;")
    println(io, "  node [fontname=\"Helvetica\", fontsize=10];")
    println(io, "  edge [fontname=\"Helvetica\", fontsize=9];")

    for nodeIndex in 1:length(graph.nodes)
        label = _dot_escape(_scenario_names(graph, nodeIndex))
        kind = graph.kinds[nodeIndex]
        if kind === :world
            println(io, "  n", nodeIndex, " [label=\"current world\\n", label,
                    "\", shape=box, style=\"bold,rounded\", penwidth=2];")
        elseif kind === :cycle
            println(io, "  n", nodeIndex, " [label=\"cycle\\n", label,
                    "\", shape=ellipse, style=dashed];")
        else
            robustnessText = isnan(graph.robustness[nodeIndex]) ? "" :
                             "\\n" * _dot_escape(_pct(graph.robustness[nodeIndex])) * " robust"
            println(io, "  n", nodeIndex, " [label=\"", label, robustnessText,
                    "\", shape=box];")
        end
    end

    for edge in graph.edges
        if isempty(edge.changes)
            println(io, "  n", edge.from, " -> n", edge.to,
                    " [label=\"(default)\", style=bold];")
        else
            label = _dot_escape(_changes_text(graph, edge; separator="\n"))
            style = length(edge.changes) > 1 ? ", style=dotted" : ""
            println(io, "  n", edge.from, " -> n", edge.to,
                    " [label=\"", replace(label, "\n" => "\\n"), "\"", style, "];")
        end
    end

    println(io, "}")
    return String(take!(io))
end

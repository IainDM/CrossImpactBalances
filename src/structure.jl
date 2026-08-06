# ══ Structure: who actually influences whom, and cutting along the seams ═══
#
# A CIB matrix holds, for every ordered pair of descriptors (i, j), a block of
# judgments: one row per variant of i, saying how that variant scores each
# variant of j. Descriptor i has a real SAY over j only when those rows
# DIFFER — if every variant of i casts identical votes on j (the all-zero
# block being the common case), then j's scores come out the same whatever i
# does, and i's state is simply inaudible to j.
#
# Draw an arrow i → j exactly where i has a say, and matrices that look like
# one big tangle often turn out to have visible structure:
#
#   ISLANDS. If the descriptors split into groups with no arrow crossing
#   between groups in either direction, each group is a self-contained CIB
#   model: under succession its descriptors read only their own group's
#   scores. Consistent scenarios of the whole model are then exactly the
#   combinations of per-group consistent ones, and — the point of this file —
#   basin sizes MULTIPLY: a start drains to the pair (X, Y) precisely when
#   its group-1 half drains to X and its group-2 half drains to Y, so the
#   count is the product of the per-group counts. Exact basin analysis of a
#   10^18 space, assembled from pieces that fit in memory (product_basins).
#
#   DIALS. A descriptor NO arrow points at (own block included) has constant
#   scores in every scenario. If those constants are all equal — nobody votes
#   on it at all — the tie rule keeps it wherever it started: it is not a
#   moving part but an external dial. Pin it to each of its variants in turn
#   (fix_descriptor) and analyse each pinned model separately: every pinned
#   analysis is a space DIVIDED by that variant count, the runs are fully
#   independent, and their results simply sit side by side. Also exact.
#
#   JUMP-AND-STAYS. If the constants are NOT all equal, the descriptor jumps
#   into its best-scoring variants on the first step and stays there. That
#   pins down where consistent scenarios can live (only at its maximizers) —
#   but it does NOT license the dial trick for basins: the starting value
#   still shoved the descriptor's neighbours around during that first step,
#   so the counts do not split. Reported for what it is, used for nothing
#   more.
#
# One tempting case is deliberately absent: group A feeding group B without
# feedback. Fixed points would split; basin COUNTS do not (where B ends up
# depends on the wobbling path A took while settling, not just on where A
# settled), so this file makes no claim there. Everything it does claim —
# island products, dial slicing — is exact, and the property tests hold it to
# the brute-force oracles.

"""
    InfluenceStructure

The influence map of a model, from [`influence_structure`](@ref). Descriptors
are identified by 1-based position in `cib.descriptors` (the display shows
names).

- `descriptors`: the names, for printing
- `activeEdges`: `activeEdges[i, j]` is true when descriptor i has a real say
  over j — some two variants of i score j's variants differently. The
  diagonal is a descriptor's say over itself (nonstandard matrices have one).
- `components`: the islands — groups of descriptors with no influence in or
  out, each sorted; singletons included. One island = no decomposition.
- `tieFrozen`: dials — descriptors nothing influences and whose constant
  scores are level, so under the tie rule they never move from where they
  start. Candidates for [`fix_descriptor`](@ref) slicing.
- `forced`: jump-and-stays — `(descriptor, maximizers)` pairs: nothing
  influences the descriptor but its constant scores have winners, so after
  one step it sits at a maximizer forever (0-based variant indices; every
  consistent scenario must hold one of them).
"""
struct InfluenceStructure
    descriptors::Vector{String}
    activeEdges::Matrix{Bool}
    components::Vector{Vector{Int}}
    tieFrozen::Vector{Int}
    forced::Vector{Tuple{Int,Vector{Int}}}
end

"""
    influence_structure(cib) -> InfluenceStructure

Compute the influence map: which descriptors actually influence which (see
the file header), the independent islands it splits into, and the dial /
jump-and-stay descriptors. Cost is one pass over the matrix — trivial even
for models whose scenario space is astronomically large, which is exactly
when it matters: the map is what decides whether [`product_basins`](@ref) or
[`fix_descriptor`](@ref) can make an impossible exact analysis possible.
"""
function influence_structure(cib::CIB)
    numberOfDescriptors = cib.numberOfDescriptors
    activeEdges = falses(numberOfDescriptors, numberOfDescriptors)
    for source in 1:numberOfDescriptors, target in 1:numberOfDescriptors
        activeEdges[source, target] = _active_edge(cib, source, target)
    end

    # Islands: connected groups of the arrows-ignoring-direction picture
    # (self-arrows don't connect anything). Plain breadth-first search.
    componentOf = zeros(Int, numberOfDescriptors)
    components = Vector{Vector{Int}}()
    queue = Int[]
    for seed in 1:numberOfDescriptors
        componentOf[seed] == 0 || continue
        push!(components, Int[])
        componentIndex = length(components)
        componentOf[seed] = componentIndex
        push!(queue, seed)
        while !isempty(queue)
            current = pop!(queue)
            push!(components[componentIndex], current)
            for other in 1:numberOfDescriptors
                other == current && continue
                componentOf[other] == 0 || continue
                if activeEdges[current, other] || activeEdges[other, current]
                    componentOf[other] = componentIndex
                    push!(queue, other)
                end
            end
        end
        sort!(components[componentIndex])
    end

    # Dials and jump-and-stays: descriptors with no incoming say (own block
    # included), classified by whether their constant scores are level.
    tieFrozen = Int[]
    forced = Tuple{Int,Vector{Int}}[]
    for descriptor in 1:numberOfDescriptors
        any(activeEdges[source, descriptor] for source in 1:numberOfDescriptors) && continue
        constantScores = _constant_block_scores(cib, descriptor)
        topScore = maximum(constantScores)
        if all(==(topScore), constantScores)
            push!(tieFrozen, descriptor)
        else
            maximizers = [variant - 1 for variant in eachindex(constantScores)
                          if constantScores[variant] == topScore]
            push!(forced, (descriptor, maximizers))
        end
    end

    return InfluenceStructure(copy(cib.descriptors), activeEdges, components,
                              tieFrozen, forced)
end

# Does descriptor `source` have a real say over `target`? True iff two of
# source's variant rows differ somewhere in target's columns.
function _active_edge(cib::CIB, source::Int, target::Int)
    cib.numberOfVariants[source] <= 1 && return false   # one variant: nothing to vary
    firstRow = cib.desc_offsets[source] + 1
    lastRow = cib.desc_offsets[source] + cib.numberOfVariants[source]
    firstColumn = cib.desc_offsets[target] + 1
    lastColumn = cib.desc_offsets[target] + cib.numberOfVariants[target]
    @inbounds for column in firstColumn:lastColumn, row in firstRow+1:lastRow
        cib.cim[row, column] != cib.cim[firstRow, column] && return true
    end
    return false
end

# The scores of `descriptor`'s variants when nothing influences it: every
# contributor's rows agree on this block, so variant-0 rows stand for all.
function _constant_block_scores(cib::CIB, descriptor::Int)
    firstColumn = cib.desc_offsets[descriptor] + 1
    lastColumn = cib.desc_offsets[descriptor] + cib.numberOfVariants[descriptor]
    scores = zeros(Int, cib.numberOfVariants[descriptor])
    for source in 1:cib.numberOfDescriptors
        representativeRow = cib.desc_offsets[source] + 1
        for (position, column) in enumerate(firstColumn:lastColumn)
            scores[position] += cib.cim[representativeRow, column]
        end
    end
    return scores
end

function Base.show(io::IO, ::MIME"text/plain", structure::InfluenceStructure)
    numberOfDescriptors = length(structure.descriptors)
    edgeCount = count(structure.activeEdges)
    println(io, "InfluenceStructure: ", numberOfDescriptors, " descriptors, ",
            edgeCount, " active influences")
    if length(structure.components) == 1
        println(io, "  Independent groups: 1 (no decomposition — every descriptor is ",
                "connected to every other, directly or indirectly)")
    else
        println(io, "  Independent groups: ", length(structure.components),
                " — exact basin analysis composes across them (product_basins):")
        for (index, component) in enumerate(structure.components)
            names = join(structure.descriptors[component], ", ")
            println(io, "    group ", index, " (", length(component),
                    length(component) == 1 ? " descriptor): " : " descriptors): ", names)
        end
    end
    if !isempty(structure.tieFrozen)
        names = join(structure.descriptors[structure.tieFrozen], ", ")
        println(io, "  Dials (never move; slice exactly with fix_descriptor): ", names)
    end
    for (descriptor, maximizers) in structure.forced
        println(io, "  Jump-and-stay: ", structure.descriptors[descriptor],
                " settles into variant(s) ", maximizers,
                " after one step (consistent scenarios must hold one)")
    end
    isempty(structure.tieFrozen) && isempty(structure.forced) &&
        println(io, "  No dial or jump-and-stay descriptors")
    print(io, "  (arrows: activeEdges[i, j] = descriptor i has a real say over j)")
end

"""
    fix_descriptor(cib, descriptor, variant) -> CIB

A new model with `descriptor` pinned to `variant`: the descriptor keeps its
place but retains only that one variant (the engine treats one-variant
descriptors as never moving), and every judgment row and column belonging to
its other variants is dropped. Both may be given by name or 0-based index, as
in [`get_impact`](@ref).

The pinned model's scenario space is the original's sliced at that setting —
smaller by the descriptor's variant count — and its dynamics are exactly the
original's restricted to that slice *when nothing gave the pinned descriptor
a reason to move* (it was a dial — see [`InfluenceStructure`](@ref) — or the
slice is simply the question being asked). Two uses:

- Conditional analysis: "hold Trade at Protectionist and show me the
  consistent scenarios / basins of everything else." Cheap to build, so
  sweeping a judgment or a setting across its range is a loop, not a
  re-parse.
- Exact decomposition: pin a dial to each of its variants in turn and run the
  pinned analyses side by side; pin everything outside an island to reduce
  the model to that island ([`split_cib`](@ref) does this).

The result is a complete, ordinary `CIB` (kernel left empty — the original's
consistent scenarios use the unpinned space's numbering) and works with every
analysis in the package.
"""
function fix_descriptor(cib::CIB, descriptor, variant)
    descriptorPosition = _desc_index(cib, descriptor)
    keptRow = _table_index(cib, descriptor, variant)   # validates the variant too
    blockStart = cib.desc_offsets[descriptorPosition]
    blockWidth = cib.numberOfVariants[descriptorPosition]

    # The rows/columns that survive: everything except the pinned descriptor's
    # block, plus the one kept variant. Slicing rows and columns by the SAME
    # index list keeps the matrix square and every other block intact.
    keptIndices = vcat(1:blockStart, keptRow, blockStart+blockWidth+1:cib.numberOfDimensions)
    newCim = cib.cim[keptIndices, keptIndices]

    newVariantCounts = copy(cib.numberOfVariants)
    newVariantCounts[descriptorPosition] = 1

    descriptorName = cib.descriptors[descriptorPosition]
    keptVariantName = cib.variants[descriptorName][keptRow - blockStart]
    newVariantNames = copy(cib.variants)               # per-descriptor lists are shared, not edited
    newVariantNames[descriptorName] = [keptVariantName]

    newOffsets = Vector{Int}(undef, cib.numberOfDescriptors)
    runningOffset = 0
    for index in 1:cib.numberOfDescriptors
        newOffsets[index] = runningOffset
        runningOffset += newVariantCounts[index]
    end

    # The transpose is REGENERATED from the sliced matrix, never sliced from
    # the old transpose — one permutedims keeps the two provably in step.
    return CIB(copy(cib.descriptors), newVariantNames, newVariantCounts,
               newCim, permutedims(newCim), size(newCim, 1),
               cib.numberOfDescriptors, Vector{Vector{Int}}(), newOffsets)
end

"""
    split_cib(cib; structure=influence_structure(cib)) -> Vector{CIB}

One model per island of the influence map: island k's model is the original
with every descriptor OUTSIDE the island pinned ([`fix_descriptor`](@ref)) to
its first variant. Pinning is exact here because no outside descriptor has
any say over the island — its rows all agree on the island's columns, so
WHICH variant stays makes no difference to the island's scores; the one kept
row carries the (variant-independent) contribution those descriptors were
always making. Descriptor positions are unchanged, so an island model's
scenarios read directly as the island's part of a full scenario.
"""
function split_cib(cib::CIB; structure::InfluenceStructure=influence_structure(cib))
    length(structure.descriptors) == cib.numberOfDescriptors || throw(ArgumentError(
        "split_cib: this InfluenceStructure describes a different model " *
        "($(length(structure.descriptors)) descriptors vs $(cib.numberOfDescriptors))"))
    return map(structure.components) do component
        reduced = cib
        for descriptor in 1:cib.numberOfDescriptors
            descriptor in component && continue
            reduced = fix_descriptor(reduced, descriptor - 1, 0)   # 0-based index form
        end
        reduced
    end
end

"""
    ComposedBasins

Exact basin analysis of a decomposable model, from [`product_basins`](@ref) —
the same information as [`find_basins`](@ref)' tuple, in `Int128` because
composed counts live at scales `Int64` cannot hold.

- `fixedPoints`, `basinSizes`: every consistent scenario of the FULL model
  (one combination per choice of island consistent scenarios), sorted by
  ascending signature, with its exact full-space basin size
- `cycleCount`: full-space starts that end in a cycle (some island cycling)
- `scenarioCount`: the full space, so `sum(basinSizes) + cycleCount == scenarioCount`
- `components`: the islands (1-based descriptor positions), and
  `componentResults`: each island's own `(fixedPoints, basinSizes, cycleCount)`
  exactly as `find_basins` returned it
"""
struct ComposedBasins
    fixedPoints::Vector{Vector{Int}}
    basinSizes::Vector{Int128}
    cycleCount::Int128
    scenarioCount::Int128
    components::Vector{Vector{Int}}
    componentResults::Vector{Tuple{Vector{Vector{Int}},Vector{Int},Int}}
end

"""
    product_basins(cib; rule=GlobalSuccession(), method=:auto,
                   structure=influence_structure(cib), cache_bytes=2^30,
                   progress=false) -> ComposedBasins

Exact basins for a model whose influence map splits into independent islands:
run [`find_basins`](@ref) on each island's model ([`split_cib`](@ref)) and
compose — island results multiply, because a start drains to a combination of
island attractors exactly when each of its island parts drains to its own
(see the file header). With islands of sizes n₁, n₂, … this prices an exact
analysis of the n₁·n₂·… -scenario model at the cost of the largest island —
the only exact route once the full space passes ~10¹³.

Composition is proven for the built-in rules, whose steps act on each island
separately; a custom rule is refused rather than guessed about. With a single
island this degenerates to plain `find_basins` (in `Int128` clothing) — the
result is then exactly as unreachable as before, which the influence map will
already have told you. `method`, `cache_bytes` and `progress` are forwarded
to each island's `find_basins`.
"""
function product_basins(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                        method::Symbol=:auto,
                        structure::InfluenceStructure=influence_structure(cib),
                        cache_bytes::Integer=1 << 30, progress::Bool=false)
    rule isa Union{GlobalSuccession,SequentialSuccession} || throw(ArgumentError(
        "product_basins: composition across islands is proven for the built-in rules " *
        "only — a custom rule need not act on each island independently. Run find_basins " *
        "on the split_cib parts yourself if your rule does."))

    islandModels = split_cib(cib; structure=structure)
    islandResults = [find_basins(model; rule=rule, method=method,
                                 cache_bytes=cache_bytes, progress=progress)
                     for model in islandModels]

    totalScenarios = scenario_count(cib)
    combinationCount = prod(Int128[length(result[1]) for result in islandResults])
    combinationCount <= 100_000 || throw(ArgumentError(
        "product_basins: the islands combine to $(combinationCount) consistent scenarios " *
        "— too many to enumerate usefully; work with componentResults via split_cib instead"))

    # Every combination of island consistent scenarios, composed positionally:
    # descriptor positions are preserved by split_cib, so each full scenario is
    # assembled by taking each island's values at its own descriptors.
    numberOfIslands = length(islandResults)
    fixedPoints = Vector{Vector{Int}}()
    basinSizes = Vector{Int128}()
    if all(result -> !isempty(result[1]), islandResults)
        choice = ones(Int, numberOfIslands)            # an odometer over island kernels
        while true
            fullScenario = Vector{Int}(undef, cib.numberOfDescriptors)
            combinedSize = Int128(1)
            for islandIndex in 1:numberOfIslands
                islandScenario = islandResults[islandIndex][1][choice[islandIndex]]
                for descriptor in structure.components[islandIndex]
                    fullScenario[descriptor] = islandScenario[descriptor]
                end
                combinedSize *= islandResults[islandIndex][2][choice[islandIndex]]
            end
            push!(fixedPoints, fullScenario)
            push!(basinSizes, combinedSize)
            position = 1
            while position <= numberOfIslands
                choice[position] += 1
                choice[position] <= length(islandResults[position][1]) && break
                choice[position] = 1
                position += 1
            end
            position > numberOfIslands && break
        end
    end

    # A full-space start converges iff every island part converges, so the
    # converging mass is the product of per-island converging masses; the rest
    # — some island cycling — is the composed cycle count.
    convergingMass = prod(Int128[sum(Int128, result[2]; init=Int128(0))
                                 for result in islandResults])
    cycleCount = totalScenarios - convergingMass

    sortOrder = sortperm(fixedPoints; by=scenario -> _signature128(cib, scenario))
    return ComposedBasins(fixedPoints[sortOrder], basinSizes[sortOrder], cycleCount,
                          totalScenarios, structure.components, islandResults)
end

function Base.show(io::IO, ::MIME"text/plain", composed::ComposedBasins)
    islandSizes = join([string(length(component)) for component in composed.components], "+")
    println(io, "ComposedBasins: ", _group_digits(composed.scenarioCount),
            " scenarios, composed exactly from ", length(composed.components),
            " independent islands (", islandSizes, " descriptors)")
    println(io, "  ", length(composed.fixedPoints), " consistent scenarios; cycles: ",
            _group_digits(composed.cycleCount), " (",
            _pct(Float64(composed.cycleCount) / Float64(composed.scenarioCount)), ")")
    order = sortperm(composed.basinSizes; rev=true)
    shown = first(order, 10)
    for index in shown
        println(io, "  ", composed.fixedPoints[index], ": ",
                _group_digits(composed.basinSizes[index]), " (",
                _pct(Float64(composed.basinSizes[index]) / Float64(composed.scenarioCount)), ")")
    end
    length(order) > length(shown) &&
        println(io, "  … and ", length(order) - length(shown), " more")
    print(io, "  sum(basinSizes) + cycleCount == scenarioCount (exact)")
end

Base.show(io::IO, composed::ComposedBasins) =
    print(io, "ComposedBasins(", length(composed.fixedPoints), " consistent scenarios, ",
          length(composed.components), " islands)")

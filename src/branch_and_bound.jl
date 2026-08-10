"""═══ How branch-and-bound works ═══════════════════════════════════════════

 THE PROBLEM
 -----------
The sweep above answers "which scenarios are consistent?" by building every single scenario and testing it. That is honest but wasteful: a model with 10 descriptors of 3 variants each has 59,049 scenarios, and one with 20 such descriptors has 3.5 billion. Usually only a handful are consistent.

THE IDEA
--------
 Instead of building whole scenarios, build them one descriptor at a time, and throw away whole families of them at once.

Picture a tree. At the top, nothing is decided. The first level picks a variant for descriptor 1, the second level for descriptor 2, and so on; each leaf at the bottom is one complete scenario. A "node" is a partly-filled-in scenario, e.g. "Economy = Boom, Policy = Green, everything else undecided".

The trick is that we can often prove a partial scenario is hopeless — that NO way of filling in the remaining descriptors could ever be consistent. When we can prove that, we discard the node and never visit anything beneath it.
Cutting a node at the point where 10 descriptors of 3 variants are still undecided throws away 3^10 = 59,049 scenarios in a single test.

HOW WE PROVE A NODE IS HOPELESS
-------------------------------
Remember what "consistent" means: for every descriptor, the variant the scenario has chosen must score at least as well as its siblings (else that descriptor would move, and the scenario would not be consistent).

So what matters is never a variant's score on its own — it is the GAP between a sibling's score and the chosen variant's score. In a partial scenario we can't know that gap exactly, because the undecided descriptors will still add to both sides. But we can bound it. The gap is:

     (what the already-decided descriptors contribute to the gap)   <- known exactly
   + (what the undecided descriptors will contribute to the gap)    <- unknown, but bounded

The first part is tracked as we descend (`prefixBalance`). For the second we precompute, once, the least each undecided descriptor could add to the gap of every (chosen, sibling) pair (`_bnb_bounds` → `sufDiff`).

The crucial point is that an undecided descriptor sits on ONE variant, and that single choice feeds the sibling and the chosen variant TOGETHER. So we do not bound the sibling's score and the chosen variant's score separately — the variant that flatters the sibling most may not be the variant that starves the chosen variant most, and pairing those two extremes imagines a world that cannot happen. Instead, for each pair we ask: which single variant of this undecided descriptor closes the gap the most? Summing that per-descriptor worst case gives the gap's guaranteed floor.

If even that guaranteed gap leaves the sibling strictly ahead (by more than the rule's margin), the sibling wins no matter how the undecided descriptors turn out. That descriptor would always want to move. Every scenario below the chosen variant is inconsistent so there's no need to look at any of them.

A WORKED EXAMPLE
----------------
Three descriptors: Economy (Boom/Bust), Policy (Green/Grey), Energy (Renewable/Fossil).
We have descended to the node

     Economy = Boom, Policy = Green, Energy = still undecided

and we are checking whether Policy is happy with Green. Between them, the descriptors decided so far contribute to Policy's two variants:

     Green: -4                      Grey: +5

so the decided part of Grey's gap over Green is +9. Energy is still undecided; its two variants contribute:

     Renewable: +1 to Green, +2 to Grey   -> widens the gap by +1
     Fossil:     0 to Green, +1 to Grey   -> widens the gap by +1

Whichever way Energy goes, the gap grows by at least min(+1, +1) = +1 — that is the precomputed `sufDiff` entry for the pair (Green, Grey). The final gap is therefore at least 9 + 1 = +10 > 0: Grey strictly beats Green in every completion, Policy abandons Green, and both scenarios under this node — (Boom, Green, Renewable) and (Boom, Green, Fossil) — are inconsistent. We discard them without ever scoring them. In a 12-descriptor model the same single test would have discarded tens of thousands.

(Bounding each column separately would have been weaker: Green ends up between -4 and -3, Grey between +6 and +7, so the guaranteed gap would only be +6 - (-3) = +9, not +10. The lost +1 comes from pairing Green's best case (Renewable) with Grey's worst case (Fossil) — two different Energy choices that can never happen at once. Here both bounds happen to prune; on tighter matrices only the pair bound does.)


WHY THE ANSWER IS STILL EXACT
-----------------------------
Two things make this a shortcut rather than an approximation:

1. The bounds are conservative. We only ever prune when the sibling stays ahead even in the completion most favourable to the chosen variant, so a scenario that could be consistent is never discarded.

2. At the bottom of the tree, when every descriptor is decided, there is nothing left undecided, so the suffix differences are all zero and the guaranteed gap IS the true gap. The prune test becomes precisely the ordinary consistency test. That is why a leaf which survives the prune IS a consistent scenario and needs no further checking.


The result is identical to the sweep's.

 WHY THERE IS A NODE BUDGET
--------------------------
Pruning only pays when the cross-impact matrix is strongly coupled, so that descriptors decisively push each other around.
If the impacts are weak or evenly balanced, the guaranteed gaps almost never clear the margin, almost nothing gets pruned, and we end up walking the whole tree — which is slower than the sweep, because the sweep has the odometer trick and visits only the leaves whereas the tree walk also visits every internal node.

testing on the sample files in the repository shows that the big, strongly-coupled models are where the win is: CIB_nested finds its 20 consistent scenarios after looking at one twenty-five-thousandth of the space.

In some cases - tiny or weakly-coupled models - this branch and bound takes more effort than it saves. So we say that it can run until it reaches a certain % of nodes and if it hasn't worked by then, it is abandoned for a full sweep using _find_kernel.

══════════════════════════════════════════════════════════════════════════
"""

"""
TECHNICAL description of _bnb_bounds calculations

    _bnb_bounds(cib) -> (sufDiff, pairOffsets)

Precompute, for every ordered pair of variant columns (chosen, rival) belonging to the same descriptor, the least the still-undecided descriptors could add to the rival's lead over the chosen variant. This is the "unknown but bounded" half of the gap described in the section header, computed once and reused by every node in the tree.

`sufDiff[p, k]` is the minimum total contribution descriptors `k..numberOfDescriptors` can make to (rival's score − chosen's score), over every possible choice of their variants. `p` identifies the pair: descriptor `i`'s pairs occupy a block starting at `pairOffsets[i]` (0-based, mirroring `desc_offsets`), and the pair (chosenLocal, rivalLocal) — both 0-based — sits at `p = pairOffsets[i] + chosenLocal * v_i + rivalLocal + 1`.

Across descriptors the minima simply add, because each undecided descriptor chooses its variant independently of the others. Within one descriptor they do NOT decompose per column: its single variant choice feeds the chosen and the rival columns together, so the minimum is taken over that shared choice — `min_a (cim[a→rival] − cim[a→chosen])`. That coupling makes this bound strictly tighter than bounding each column's score separately whenever different variants attain the two per-column extremes.

Column `numberOfDescriptors + 1` is left as zero — "no descriptors left to decide, so they can add nothing". That is what makes the prune test exact at the bottom of the tree. So once descriptors `1..k` are assigned, `sufDiff[p, k+1]` is the guaranteed floor on what the still-undecided descriptors add to pair `p`'s gap.

The table is stored pairs-major, `(number of pairs) × (numberOfDescriptors + 1)`: Julia arrays are column-major, so the prune's scan over one descriptor's rivals at a fixed depth reads consecutive memory.

A self-pair (rival == chosen) has difference identically zero, so its entries are all zero. It is kept so each descriptor's rival scan is one contiguous block, and so that a (hypothetical) negative margin reproduces the sweep's semantics exactly — the sweep's rival loop includes the chosen column too.
"""
function _bnb_bounds(cib::CIB)
    numberOfDescriptors = cib.numberOfDescriptors
    cimTranspose = cib.cim_t

    pairOffsets = Vector{Int}(undef, numberOfDescriptors)
    numberOfPairs = 0
    for descriptorIndex in 1:numberOfDescriptors
        pairOffsets[descriptorIndex] = numberOfPairs
        numberOfPairs += cib.numberOfVariants[descriptorIndex]^2
    end

    sufDiff = zeros(Int, numberOfPairs, numberOfDescriptors + 1)

    # Build the bounds back-to-front: the floor from descriptor k onwards is
    # its own worst single-variant gap contribution plus whatever descriptors
    # k+1..end can add. `ndesc:-1:1` is a range counting down (like
    # range(n, 0, -1) in Python).
    @inbounds for sourceDescriptor in numberOfDescriptors:-1:1
        sourceOffset = cib.desc_offsets[sourceDescriptor]
        sourceVariantCount = cib.numberOfVariants[sourceDescriptor]
        pairIndex = 0
        for descriptorIndex in 1:numberOfDescriptors
            columnOffset = cib.desc_offsets[descriptorIndex]
            variantCount = cib.numberOfVariants[descriptorIndex]
            for chosenLocal in 0:variantCount-1
                chosenColumn = columnOffset + chosenLocal + 1
                for rivalLocal in 0:variantCount-1
                    rivalColumn = columnOffset + rivalLocal + 1
                    pairIndex += 1
                    minDiff = typemax(Int)
                    for sourceVariant in 0:sourceVariantCount-1
                        # One variant of the source descriptor feeds BOTH
                        # columns, so take the extreme of the difference —
                        # not the difference of the extremes.
                        sourceRow = sourceOffset + sourceVariant + 1
                        difference = Int(cimTranspose[rivalColumn, sourceRow]) -
                                     Int(cimTranspose[chosenColumn, sourceRow])
                        minDiff = ifelse(difference < minDiff, difference, minDiff)
                    end
                    sufDiff[pairIndex, sourceDescriptor] =
                        sufDiff[pairIndex, sourceDescriptor + 1] + minDiff
                end
            end
        end
    end
    return sufDiff, pairOffsets
end

"""
Everything one task needs while walking its part of the tree, bundled into a struct so the recursion passes one argument instead of a dozen.

`partialScenario` is the node we are currently at — the variants decided so far. `prefixBalance` is what those decided descriptors contribute to every variant's score: the "known exactly" half of the gap in the section header, kept up to date as we descend rather than recomputed at each node.

The three fields after `foundFixedPoints` implement the node budget. `nodesVisited` is this task's private counter (a `Ref` is a mutable single-value box, needed because the struct itself is immutable), while `totalNodes` and `abortFlag` are Atomics shared by every task — the Julia equivalent of C#'s Interlocked operations.
"""
struct _BnBState
    cim_t::Matrix{Int}
    sufdiff::Matrix{Int}
    pairOffsets::Vector{Int}
    nvariants::Vector{Int}
    offsets::Vector{Int}
    ndesc::Int
    ndim::Int
    partialScenario::Vector{Int}
    prefixBalance::Vector{Int}
    foundFixedPoints::Vector{Vector{Int}}
    nodesVisited::Base.RefValue{Int}
    totalNodes::Threads.Atomic{Int}
    abortFlag::Threads.Atomic{Bool}
    nodeBudget::Int
    margin::Int
end

"""
    _bnb_pruned(state, assignedCount) -> Bool

Can this whole branch of the tree be thrown away? This is the test worked through in the section header: with descriptors `1..assignedCount` decided, return true if some decided descriptor has a rival variant whose guaranteed gap over the chosen variant — the decided part plus the precomputed floor on what the undecided descriptors add — exceeds the rule's `margin`. Such a rival wins in *every* possible completion, so no scenario below this node can be consistent and the caller skips the entire subtree.

Only *decided* descriptors are checked, because an undecided one has no chosen variant to be unhappy with yet.

The rival scan includes the chosen variant itself (its pair entries are all zero, so for `margin ≥ 0` it can never trigger); keeping it makes the scan one contiguous block and preserves exact sweep parity for negative margins.

Ties — and gaps within the margin — never prune, which matches the convention everywhere else that a descriptor keeps its current variant unless something strictly beats it. Once every descriptor is assigned the suffix differences are zero, so this condition becomes exactly the margin consistency test on the completed scenario (`margin = 0` recovers plain global consistency).
"""
function _bnb_pruned(state::_BnBState, assignedCount::Int)
    prefixBalance = state.prefixBalance
    sufDiff = state.sufdiff
    suffixColumn = assignedCount + 1   # bounds column for the still-undecided descriptors
    @inbounds for descriptorIndex in 1:assignedCount
        offset = state.offsets[descriptorIndex]
        variantCount = state.nvariants[descriptorIndex]
        chosenLocal = state.partialScenario[descriptorIndex]
        # What a rival must overcome: the chosen variant's decided score plus
        # the margin. In the worked example this is Green's -4.
        threshold = prefixBalance[offset + chosenLocal + 1] + state.margin
        pairBase = state.pairOffsets[descriptorIndex] + chosenLocal * variantCount
        for rivalLocal in 1:variantCount
            # The rival's decided score plus the guaranteed floor on what the
            # undecided descriptors add to its lead — in the example, Grey's
            # +5 + 1 = +6 against Green's -4. One number per PAIR, because
            # each undecided descriptor feeds both columns with the same
            # variant choice.
            if prefixBalance[offset + rivalLocal] + sufDiff[pairBase + rivalLocal, suffixColumn] > threshold
                return true
            end
        end
    end
    return false
end

"""
Count one visited node against the budget (see "why there is a node budget" in the section header). Every task keeps its own private tally and only adds
# it to the shared total every 256 nodes, so the threads rarely have to touch shared state — the exact count doesn't matter, only whether we are roughly
# over budget. Returns false once the budget has been blown, which tells the whole search to give up and let the caller fall back to the sweep.
"""
function _bnb_charge!(state::_BnBState)
    # `state.nodesVisited[]` — the [] reads/writes the value inside a Ref box.
    state.nodesVisited[] += 1
    if state.nodesVisited[] >= 256
        Threads.atomic_add!(state.totalNodes, state.nodesVisited[])
        state.nodesVisited[] = 0
        if state.totalNodes[] > state.nodeBudget || state.abortFlag[]
            state.abortFlag[] = true
            return false
        end
    end
    return true
end

"""
Visit one level of the tree: try each variant of descriptor `depth` in turn, and for each one either prune it or recurse into the descriptors below it. This is the depth-first walk itself — "decide one more descriptor, see if the result is already hopeless, and if not keep going".

Only one scenario buffer exists per task, so each variant is pushed onto the running prefix balance before descending and popped off afterwards, leaving the state exactly as it was found. That is why there is no per-node allocation despite the recursion.

Returns false if the node budget ran out and the search is abandoning.
"""
function _bnb_node!(state::_BnBState, depth::Int)
    variantCount = state.nvariants[depth]
    offset = state.offsets[depth]
    prefixBalance = state.prefixBalance
    cimTranspose = state.cim_t
    numberOfDimensions = state.ndim
    @inbounds for variantIndex in 0:variantCount-1
        variantColumn = offset + variantIndex + 1
        # Push this variant's row onto the running prefix balance...
        @simd for targetVariant in 1:numberOfDimensions
            prefixBalance[targetVariant] += cimTranspose[targetVariant, variantColumn]
        end
        state.partialScenario[depth] = variantIndex
        keepGoing = _bnb_charge!(state)
        if keepGoing && !_bnb_pruned(state, depth)
            if depth == state.ndesc
                # All descriptors assigned and the prune (now exact) passed:
                # this scenario is consistent. Store a copy — the buffer is reused.
                push!(state.foundFixedPoints, copy(state.partialScenario))
            else
                keepGoing = _bnb_node!(state, depth + 1)
            end
        end
        # ...and pop it again before trying the next variant.
        @simd for targetVariant in 1:numberOfDimensions
            prefixBalance[targetVariant] -= cimTranspose[targetVariant, variantColumn]
        end
        keepGoing || return false
    end
    return true
end

"""
    _bnb_fixed_points(cib, sufDiff, pairOffsets; node_budget, margin)
        -> (Union{Nothing, Vector{Vector{Int}}}, nodesVisited)

Run the branch-and-bound search described in the section header, across all threads.

The work is split by chopping the tree near the top: every combination of variants for the first few descriptors becomes one task, which then explores its own subtree alone. Enough descriptors are taken to make comfortably more tasks than threads, so no thread sits idle waiting for a slow branch. Tasks share nothing but the node counter, and each collects its own results.

Returns a tuple. The first element is the complete kernel sorted by ascending signature — ordered by comparing variant choices from the last descriptor down, so the order is right even for a model whose signatures overrun `Int64` — or `nothing` if the search gave up because it blew the node budget, which is the caller's signal to fall back to the sweep. The second is how many tree nodes were visited, i.e. how many partial scenarios were expanded; comparing that with the total number of scenarios shows how much the pruning actually saved.
"""
function _bnb_fixed_points(cib::CIB, sufDiff::Matrix{Int}, pairOffsets::Vector{Int};
                           node_budget::Int, margin::Int=0)
    numberOfDescriptors = cib.numberOfDescriptors
    variantCounts = cib.numberOfVariants
    descriptorOffsets = cib.desc_offsets
    numberOfDimensions = cib.numberOfDimensions
    cimTranspose = cib.cim_t

    # Choose how many leading descriptors to decide up front. Every combination
    # of their variants becomes one task, so keep taking descriptors until
    # there are comfortably more tasks than threads. E.g. with 4 threads and
    # 3-variant descriptors, two descriptors give 9 tasks, three give 27.
    prefixDepth = 0
    numberOfPrefixes = 1
    while prefixDepth < numberOfDescriptors && numberOfPrefixes < 4 * Threads.nthreads()
        prefixDepth += 1
        numberOfPrefixes *= variantCounts[prefixDepth]
    end

    totalNodes = Threads.Atomic{Int}(0)
    abortFlag = Threads.Atomic{Bool}(false)
    taskOutputs = [Vector{Vector{Int}}() for _ in 1:numberOfPrefixes]

    @sync for prefixId in 0:numberOfPrefixes-1
        taskOutput = taskOutputs[prefixId + 1]
        Threads.@spawn begin
            state = _BnBState(cimTranspose, sufDiff, pairOffsets, variantCounts, descriptorOffsets,
                              numberOfDescriptors, numberOfDimensions,
                              zeros(Int, numberOfDescriptors), zeros(Int, numberOfDimensions),
                              taskOutput, Ref(0), totalNodes, abortFlag, node_budget, margin)
            # Decode this task's prefixId into variant choices for the first
            # prefixDepth descriptors (same mixed-radix arithmetic as
            # inv_signature), checking the prune at every level so a subtree
            # already dead partway through the prefix is skipped without
            # descending into it.
            remainder = prefixId
            alive = true
            @inbounds for descriptorIndex in 1:prefixDepth
                state.partialScenario[descriptorIndex] = remainder % variantCounts[descriptorIndex]
                remainder = remainder ÷ variantCounts[descriptorIndex]
                variantColumn = descriptorOffsets[descriptorIndex] +
                                state.partialScenario[descriptorIndex] + 1
                @simd for targetVariant in 1:numberOfDimensions
                    state.prefixBalance[targetVariant] += cimTranspose[targetVariant, variantColumn]
                end
                alive = _bnb_charge!(state) && !_bnb_pruned(state, descriptorIndex)
                alive || break
            end
            if alive
                if prefixDepth == numberOfDescriptors
                    # The prefix IS a full scenario and it survived the
                    # (exact, at full depth) prune: it is consistent.
                    push!(taskOutput, copy(state.partialScenario))
                else
                    _bnb_node!(state, prefixDepth + 1)
                end
            end
            Threads.atomic_add!(totalNodes, state.nodesVisited[])   # flush the residual private count
        end
    end

    # `abortFlag[]` reads the atomic's value.
    abortFlag[] && return (nothing, totalNodes[])
    kernel = Vector{Vector{Int}}()
    for taskOutput in taskOutputs
        append!(kernel, taskOutput)
    end
    # Tasks finish in nondeterministic order; sort so the result is stable.
    # `_sig_less` (walk.jl) compares the variant choices from the LAST
    # descriptor downwards, which IS ascending signature order — a signature is
    # a mixed-radix number whose most significant digit is the last descriptor —
    # while computing no signature at all. That matters twice. It is cheaper
    # than the old key, which `sort!` recomputed at every comparison rather than
    # once per element. And it is the only version that is CORRECT past
    # typemax(Int64): there `signature` wraps, so on an 80-binary-descriptor
    # model 18,432 distinct scenarios collapse onto 9,024 distinct keys, and the
    # all-ones scenario keys as a negative number and sorts FIRST.
    sort!(kernel; lt = _sig_less)
    return (kernel, totalNodes[])
end

# ══ Which descriptor to decide first ═══════════════════════════════════════
#
# The tree above is built in whatever order the descriptors happen to appear in
# the file, and that order is not neutral — it decides how soon a node can be
# proved hopeless.
#
# A node is pruned by comparing gaps that only the ALREADY-DECIDED descriptors
# contribute to. Decide five descriptors that barely touch each other and
# `prefixBalance` is still near zero: nothing to prune with, so the search
# grinds down to the leaves. Decide five that push hard on one another and the
# gaps are already wide at depth five, killing whole subtrees at the top of the
# tree where they are biggest.
#
# So: put the strongly-coupled descriptors first. Start from the descriptor
# with the most impact weight in total, then keep taking whichever remaining
# descriptor is most tied to the ones already placed. This is a cheap greedy
# heuristic, not an optimum, and it cannot affect the ANSWER — it is undone
# before the kernel is returned. It only affects how long the answer takes.
#
# On Weimer-Jehle's corpus the difference is not marginal. D60 (2.2e23
# scenarios) went from unfinished after 1.5e9 nodes to 385,394 nodes; B80
# (1.2e24) from 254 million nodes to 215,530. The repository's own fixtures all
# improve too, by between 1.06x and 4.7x — the heuristic has not been observed
# to make any model worse.

"""
    _bnb_descriptor_order(cib) -> Vector{Int}

The order to assign descriptors in, most-coupled first. `order[k]` is the original index of the descriptor branch-and-bound should decide at depth `k`. See the section header above for why the order matters.

Deterministic: ties go to the lowest original index, so the result never depends on thread count, and a model whose matrix is entirely zero comes back in file order.
"""
function _bnb_descriptor_order(cib::CIB)
    numberOfDescriptors = cib.numberOfDescriptors
    numberOfDescriptors <= 2 && return collect(1:numberOfDescriptors)

    offsets = cib.desc_offsets
    variantCounts = cib.numberOfVariants

    # coupling[i, j] = how much impact weight descriptor i's variants place on
    # descriptor j's, summed over the block and taken as absolute value: a
    # strong push towards a variant constrains j exactly as much as a strong
    # push away from it. The diagonal stays zero — a descriptor's influence on
    # itself tells us nothing about which to decide first.
    coupling = zeros(Int, numberOfDescriptors, numberOfDescriptors)
    @inbounds for sourceDescriptor in 1:numberOfDescriptors,
                  targetDescriptor in 1:numberOfDescriptors
        sourceDescriptor == targetDescriptor && continue
        sourceRows = (offsets[sourceDescriptor] + 1):(offsets[sourceDescriptor] + variantCounts[sourceDescriptor])
        targetCols = (offsets[targetDescriptor] + 1):(offsets[targetDescriptor] + variantCounts[targetDescriptor])
        total = 0
        for row in sourceRows, col in targetCols
            total += abs(cib.cim[row, col])
        end
        coupling[sourceDescriptor, targetDescriptor] = total
    end

    # `tiesToPlaced[i]` is descriptor i's total coupling to everything chosen so
    # far, in both directions, kept up to date as each choice is made rather
    # than recomputed. Before anything is placed it seeds with total coupling to
    # the whole model, which picks the busiest descriptor to start from.
    tiesToPlaced = [sum(@view coupling[i, :]) + sum(@view coupling[:, i])
                    for i in 1:numberOfDescriptors]
    alreadyPlaced = falses(numberOfDescriptors)
    order = Vector{Int}(undef, numberOfDescriptors)

    for depth in 1:numberOfDescriptors
        # Strict `>` keeps the lowest index among equals, which is what makes
        # this deterministic.
        best = 0
        bestScore = -1
        @inbounds for i in 1:numberOfDescriptors
            alreadyPlaced[i] && continue
            if tiesToPlaced[i] > bestScore
                bestScore = tiesToPlaced[i]
                best = i
            end
        end
        order[depth] = best
        alreadyPlaced[best] = true

        # The seed's score meant "coupling to everything"; from here on the
        # score means "coupling to the placed set", so the first pick resets it
        # rather than adding to it.
        @inbounds for i in 1:numberOfDescriptors
            alreadyPlaced[i] && continue
            contribution = coupling[i, best] + coupling[best, i]
            tiesToPlaced[i] = depth == 1 ? contribution : tiesToPlaced[i] + contribution
        end
    end
    return order
end

"""
    _permute_descriptors(cib, order) -> CIB

The same model with its descriptors reordered, `order` giving the original index of each new position. The variant blocks of the CIM move with their descriptors, so the model is genuinely identical — only the numbering changes.

`variants` is shared, not rebuilt: it is keyed by descriptor name, so descriptor order does not enter into it. The kernel is left empty; this model exists only to be searched.
"""
function _permute_descriptors(cib::CIB, order::Vector{Int})
    variantCounts = cib.numberOfVariants[order]

    # Where each descriptor's variant block lands in the new numbering.
    newOffsets = Vector{Int}(undef, cib.numberOfDescriptors)
    runningOffset = 0
    for descriptorIndex in 1:cib.numberOfDescriptors
        newOffsets[descriptorIndex] = runningOffset
        runningOffset += variantCounts[descriptorIndex]
    end

    # The row/column permutation the CIM needs: every variant of the new first
    # descriptor, then every variant of the new second, and so on.
    variantOrder = Vector{Int}(undef, cib.numberOfDimensions)
    next = 1
    for originalDescriptor in order
        blockStart = cib.desc_offsets[originalDescriptor] + 1
        blockEnd = cib.desc_offsets[originalDescriptor] + cib.numberOfVariants[originalDescriptor]
        for column in blockStart:blockEnd
            variantOrder[next] = column
            next += 1
        end
    end

    permutedCim = cib.cim[variantOrder, variantOrder]
    return CIB(cib.descriptors[order], cib.variants, variantCounts,
               permutedCim, permutedims(permutedCim),
               cib.numberOfDimensions, cib.numberOfDescriptors,
               Vector{Vector{Int}}(), newOffsets)
end

"""
    _bnb_search(cib; node_budget, margin) -> (Union{Nothing,Vector{Vector{Int}}}, nodesVisited)

Branch-and-bound end to end: pick a descriptor order, search in it, and translate the answer back into the caller's own descriptor numbering.

This is the whole branch-and-bound entry point — [`_bnb_descriptor_order`](@ref), [`_permute_descriptors`](@ref), `_bnb_bounds` and `_bnb_fixed_points` are its parts, and callers should not have to assemble them. `nothing` means the node budget was spent without finishing, exactly as in `_bnb_fixed_points`.
"""
function _bnb_search(cib::CIB; node_budget::Int, margin::Int=0)
    order = _bnb_descriptor_order(cib)
    searchModel = _permute_descriptors(cib, order)

    sufDiff, pairOffsets = _bnb_bounds(searchModel)
    kernel, nodesVisited = _bnb_fixed_points(searchModel, sufDiff, pairOffsets;
                                             node_budget=node_budget, margin=margin)
    isnothing(kernel) && return (nothing, nodesVisited)

    # Back to the caller's numbering: the variant at new position k belongs to
    # original descriptor order[k]. The sort inside _bnb_fixed_points ordered
    # these by the SEARCH order's signatures, which is a different order, so
    # they have to be sorted again once renumbered.
    original = Vector{Vector{Int}}(undef, length(kernel))
    for (index, searchScenario) in pairs(kernel)
        scenario = Vector{Int}(undef, cib.numberOfDescriptors)
        @inbounds for depth in 1:cib.numberOfDescriptors
            scenario[order[depth]] = searchScenario[depth]
        end
        original[index] = scenario
    end
    sort!(original; lt = _sig_less)
    return (original, nodesVisited)
end

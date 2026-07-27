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

In a partial scenario we can't know a variant's final score, because the undecided descriptors will still add to it. But we can put limits on its score. A variant's final score is:

     (what the already-decided descriptors contribute)   <- known exactly
   + (what the undecided descriptors will contribute)    <- unknown, but bounded

The first part is tracked as we descend (`prefixBalance`). For the second we precompute, once, the smallest and largest total the undecided descriptors could possibly contribute (`_bnb_bounds` → `suffixMin` / `suffixMax`).

 Now compare a chosen variant with one of its siblings:

   * chosen variant's BEST case  = its prefix + the most the rest could add
   * sibling's WORST case        = its prefix + the least the rest could add

If the sibling's worst case still beats the chosen variant's best case, then the sibling wins no matter how the undecided descriptors turn out. That descriptor would always want to move. Every scenario below the chosen variant is inconsistent so there's no need to look at any of them. 

A WORKED EXAMPLE
----------------
Three descriptors: Economy (Boom/Bust), Policy (Green/Grey), Energy (Renewable/Fossil).
We have descended to the node

     Economy = Boom, Policy = Green, Energy = still undecided

and we are checking whether Policy is happy with Green. Between them, the descriptors decided so far contribute to Policy's two variants:

     Green: -4                      Grey: +5

Energy is undecided; from the precomputed bounds, it will contribute:

     to Green: between  0 and +1    to Grey: between +1 and +2

 So whatever Energy does:

     Green ends up between -4 and -3   -> its BEST possible score is  -3
     Grey  ends up between +6 and +7   -> its WORST possible score is +6

Grey's worst (+6) still beats Green's best (-3). So in every completion, Policy abandons Green for Grey. Both scenarios under this node — (Boom, Green, Renewable) and (Boom, Green, Fossil) — are inconsistent, and we discard them without ever scoring them. In a 12-descriptor model the same single test would have discarded tens of thousands.


WHY THE ANSWER IS STILL EXACT
-----------------------------
Two things make this a shortcut rather than an approximation:

1. The bounds are conservative. We only ever prune when the sibling wins even in the most favourable case for the chosen variant, so a scenario that could be consistent is never discarded.

2. At the bottom of the tree, when every descriptor is decided, there is nothing left undecided, so the suffix bounds are all zero and both limits collapse to the true score. The prune test becomes precisely the ordinary consistency test. That is why a leaf which survives the prune IS a consistent scenario and needs no further checking.


The result is identical to the sweep's.

 WHY THERE IS A NODE BUDGET
--------------------------
Pruning only pays when the cross-impact matrix is strongly coupled, so that descriptors decisively push each other around.
If the impacts are weak or evenly balanced, the brackets overlap almost everywhere, almost nothing gets pruned, and we end up walking the whole tree — which is slower than the sweep, because the sweep has the odometer trick and visits only the leaves whereas the tree walk also visits every internal node.

testing on the sample files in the repository shows that the big, strongly-coupled models are where the win is: CIB_nested finds its 20 consistent scenarios after looking at one ten-thousandth of the space.

In some cases - tiny or weakly-coupled models - this branch and bound takes more effort than it saves. So we say that it can run until it reaches a certain % of nodes and if it hasn't worked by then, it is abandoned for a full sweep using _find_kernel.

══════════════════════════════════════════════════════════════════════════
"""

"""
TECHNICAL description of _bnb_bounds calculations

    _bnb_bounds(cib) -> (suffixMin, suffixMax)

Precompute, for every variant column, how much the descriptors from `k` onwards could add to that column's score — at least (`suffixMin`) and at most (`suffixMax`). This is the "unknown but bounded" half of the bracket described in the section header, and it is computed once and reused by every node in the tree.

`suffixMin[k, c]` / `suffixMax[k, c]` is the minimum / maximum total contribution descriptors `k..numberOfDescriptors` can make to the impact score of variant column `c`, over every possible choice of their variants. Each descriptor is free to choose independently, so the extremes just add up: the minimum is the sum of each descriptor's own smallest entry in column `c`.

Row `numberOfDescriptors + 1` is left as zero — "no descriptors left to decide, so they can add nothing". That is what makes the prune test exact at the bottom of the tree. So once descriptors `1..k` are assigned, `suffixMin[k+1, c] .. suffixMax[k+1, c]` brackets what the still-undecided descriptors can add to column `c`.
"""
function _bnb_bounds(cib::CIB)
    numberOfDescriptors, numberOfDimensions = cib.numberOfDescriptors, cib.numberOfDimensions
    cimTranspose = cib.cim_t
    suffixMin = zeros(Int, numberOfDescriptors + 1, numberOfDimensions)
    suffixMax = zeros(Int, numberOfDescriptors + 1, numberOfDimensions)

    # Build the bounds back-to-front: descriptor k's bounds are its own best/worst row entry plus whatever descriptors k+1..end can add.
    # `ndesc:-1:1` is a range counting down (like range(n, 0, -1) in Python).
    @inbounds for descriptorIndex in numberOfDescriptors:-1:1
        offset = cib.desc_offsets[descriptorIndex]
        variantCount = cib.numberOfVariants[descriptorIndex]
        for targetVariant in 1:numberOfDimensions
            minValue = typemax(Int)
            maxValue = typemin(Int)
            for variantIndex in 0:variantCount-1
                # = cim[row of this variant, targetVariant]
                value = Int(cimTranspose[targetVariant, offset + variantIndex + 1])
                minValue = ifelse(value < minValue, value, minValue)
                maxValue = ifelse(value > maxValue, value, maxValue)
            end
            suffixMin[descriptorIndex, targetVariant] = suffixMin[descriptorIndex + 1, targetVariant] + minValue
            suffixMax[descriptorIndex, targetVariant] = suffixMax[descriptorIndex + 1, targetVariant] + maxValue
        end
    end
    return suffixMin, suffixMax
end

"""
Everything one task needs while walking its part of the tree, bundled into a struct so the recursion passes one argument instead of a dozen.

`partialScenario` is the node we are currently at — the variants decided so far. `prefixBalance` is what those decided descriptors contribute to every variant's score: the "known exactly" half of the bracket in the section header, kept up to date as we descend rather than recomputed at each node.

The three fields after `foundFixedPoints` implement the node budget. `nodesVisited` is this task's private counter (a `Ref` is a mutable single-value box, needed because the struct itself is immutable), while `totalNodes` and `abortFlag` are Atomics shared by every task — the Julia equivalent of C#'s Interlocked operations.
"""
struct _BnBState
    cim_t::Matrix{Int}
    sufmin::Matrix{Int}
    sufmax::Matrix{Int}
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

Can this whole branch of the tree be thrown away? This is the test worked through in the section header: with descriptors `1..assignedCount` decided, return true if some decided descriptor has a rival variant that wins in *every* possible completion — its worst case still beats the chosen variant's best case (by more than the rule's `margin`). If so, no scenario below this node can be consistent, so the caller skips the entire subtree.

Only *decided* descriptors are checked, because an undecided one has no chosen variant to be unhappy with yet.

Ties — and gaps within the margin — never prune, which matches the convention everywhere else that a descriptor keeps its current variant unless something strictly beats it. Once every descriptor is assigned the suffix bounds are zero, so this condition becomes exactly the margin consistency test on the completed scenario (`margin = 0` recovers plain global consistency).
"""
function _bnb_pruned(state::_BnBState, assignedCount::Int)
    prefixBalance = state.prefixBalance
    suffixMin = state.sufmin
    suffixMax = state.sufmax
    suffixRow = assignedCount + 1   # bounds row for the still-undecided descriptors
    @inbounds for descriptorIndex in 1:assignedCount
        offset = state.offsets[descriptorIndex]
        chosenColumn = offset + state.partialScenario[descriptorIndex] + 1
        # Best case for the variant this branch chose: what the decided descriptors already give it, plus the most the undecided ones could add. A rival has to clear this (and the margin) to win outright. 
        # In the worked example this is Green's -4 + 1 = -3.
        bestCaseChosen = prefixBalance[chosenColumn] + suffixMax[suffixRow, chosenColumn] + state.margin
        for column in offset+1:offset+state.nvariants[descriptorIndex]
            # Worst case for a rival variant: its decided contribution plus the least the undecided descriptors could add — in the example, Grey's +5 + 1 = +6. If even that beats the chosen variant's best case, the rival wins however the rest turns out, so nothing below this node can be consistent.
            if prefixBalance[column] + suffixMin[suffixRow, column] > bestCaseChosen
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
    _bnb_fixed_points(cib, suffixMin, suffixMax; node_budget, margin)
        -> (Union{Nothing, Vector{Vector{Int}}}, nodesVisited)

Run the branch-and-bound search described in the section header, across all threads.

The work is split by chopping the tree near the top: every combination of variants for the first few descriptors becomes one task, which then explores its own subtree alone. Enough descriptors are taken to make comfortably more tasks than threads, so no thread sits idle waiting for a slow branch. Tasks share nothing but the node counter, and each collects its own results.

Returns a tuple. The first element is the complete kernel sorted by ascending signature — or `nothing` if the search gave up because it blew the node budget, which is the caller's signal to fall back to the sweep. The second is how many tree nodes were visited, i.e. how many partial scenarios were expanded; comparing that with the total number of scenarios shows how much the pruning actually saved.
"""
function _bnb_fixed_points(cib::CIB, sufmin::Matrix{Int}, sufmax::Matrix{Int};
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
            state = _BnBState(cimTranspose, sufmin, sufmax, variantCounts, descriptorOffsets,
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
    sort!(kernel; by = scenario -> signature(cib, scenario))
    return (kernel, totalNodes[])
end

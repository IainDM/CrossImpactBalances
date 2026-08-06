""" ══ Basins of attraction: where does each scenario end up? ═════════════════

find_consistent asks which scenarios are stable. This asks a bigger question: start from EVERY scenario in the model, keep applying succession, and see where you land. The number of starting points that drain into a given consistent scenario is its basin — a measure of how much of the possibility space that future attracts.

A NOTE ON WORDING. This file says "fixed point" more than the rest of the codebase, and deliberately. Here we are treating succession as a dynamical system, and a chain can end in one of two ways: at a scenario that maps to itself (a fixed point — which is exactly what the rest of the code calls a consistent scenario), or in a repeating loop of two or more scenarios that chase each other forever (a cycle). Cycles are not consistent scenarios and are not anybody's basin, so they are counted separately.
"""

"""
    find_basins(cib; rule=GlobalSuccession()) -> (fixed_points, basin_sizes, cycle_count)

Exhaustive basin-of-attraction analysis under `rule` (default [`GlobalSuccession`](@ref)). Follows the succession chain from every scenario in the space, counting how many starting points converge to each consistent scenario.

Returns three things: the consistent scenarios found, how many scenarios drain into each one, and how many scenarios instead end up in a repeating cycle and so belong to no basin. (A "fixed point" here is the same thing as a consistent scenario — see the note above.)

This runs in two phases: a threaded sweep fills a flat successor table (every scenario's succession step), then a resolution pass walks the table with path compression so each scenario is resolved exactly once. Memory is two flat tables — the successor table at 4 bytes per scenario when every signature fits `Int32`, 8 otherwise, plus a 4-byte-per-scenario label table — independent of the thread count. A model too big for that memory raises an informative error rather than crashing; see `method` below for the alternatives.
For example, the first phase works out that scenario A maps to B, B maps to C, and C maps to D which is a consistent scenario. Therefore A, B, C and D are all in D's basin, and the second phase works this out.

The default [`GlobalSuccession`](@ref) rule takes a fast path: the successor table is filled by a mixed-radix odometer that maintains the impact balance incrementally (one row delta per scenario, no allocation) instead of calling [`succession_step`](@ref) per scenario. Any other rule uses the generic per-scenario path — identical output, just slower.
"""
function find_basins(cib::CIB; rule::SuccessionRule=GlobalSuccession())
    return _find_basins(rule, cib)
end

#region "Fast path"
# Fast path, selected by dispatch when the rule is exactly GlobalSuccession. This code is fine tuned to calculate that rule quickly
function _find_basins(::GlobalSuccession, cib::CIB)
    numberOfScenarios = max_signature(cib) + 1
    scoreType = _score_type(cib)
    # Store table entries as Int32 when every signature fits — half the
    # memory, and memory traffic is what bounds this analysis.
    if numberOfScenarios <= Int(typemax(Int32)) - 1
        return _fast_basins(cib, Int32, Matrix{scoreType}(cib.cim_t))
    end
    return _fast_basins(cib, Int64, Matrix{scoreType}(cib.cim_t))
end

# `::Type{SignatureInt}` receives a TYPE as an argument value (Int32 or Int64 above) — passing types around as ordinary values is normal in Julia.
function _fast_basins(cib::CIB, ::Type{SignatureInt},
                      cimTranspose::Matrix{ScoreInt}) where {SignatureInt<:Union{Int32,Int64}, ScoreInt<:Signed}
    numberOfScenarios = max_signature(cib) + 1
    successorTable = Vector{SignatureInt}(undef, numberOfScenarios)
    _successor_table!(successorTable, cib, cimTranspose)
    fixedPointSignatures, basinSizes, cycleCount = _resolve_and_tally(successorTable, numberOfScenarios) #grab the basin details
    fixedPoints = [inv_signature(cib, sig) for sig in fixedPointSignatures] #reverse out the fixed points in terms of descriptors and variants rather than index
    return (fixedPoints, basinSizes, cycleCount)
end
#endregion

#region "Generic path"
# Generic path: works for any rule, because it only ever calls the rule's own succession_step. Allocates a few vectors per scenario, so it is much slower than the odometer path — a correctness baseline or for use when there is no alternative
function _find_basins(rule::SuccessionRule, cib::CIB)
    numberOfScenarios = max_signature(cib) + 1
    if numberOfScenarios <= Int(typemax(Int32)) - 1
        return _generic_basins(rule, cib, Int32)
    end
    return _generic_basins(rule, cib, Int64)
end

function _generic_basins(rule::SuccessionRule, cib::CIB,
                         ::Type{SignatureInt}) where {SignatureInt<:Union{Int32,Int64}}
    numberOfScenarios = max_signature(cib) + 1
    successorTable = Vector{SignatureInt}(undef, numberOfScenarios)
    numberOfThreads = Threads.nthreads()
    chunkSize = cld(numberOfScenarios, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * chunkSize
        lastSignature = min(threadIndex * chunkSize, numberOfScenarios) - 1
        firstSignature > lastSignature && continue
        # Each task writes a disjoint range of the table, so there is no data race even though they share the array.
        Threads.@spawn for currentSignature in firstSignature:lastSignature
            scenario = inv_signature(cib, currentSignature)
            successor = succession_step(rule, cib, scenario)
            @inbounds successorTable[currentSignature + 1] =
                SignatureInt(signature(cib, successor))
        end
    end
    fixedPointSignatures, basinSizes, cycleCount = _resolve_and_tally(successorTable, numberOfScenarios)
    fixedPoints = [inv_signature(cib, sig) for sig in fixedPointSignatures]
    return (fixedPoints, basinSizes, cycleCount)
end
#endregion

#region "Successor table"

"""
══ Phase 1: building the successor table ══════════════════════════════════

The successor table is a single flat array holding, for every scenario in the model, the number of the scenario it steps to:

     successorTable[sig + 1] = the signature of sig's successor

(The `+ 1` throughout is only Julia's 1-based indexing meeting our 0-based scenario numbering. Scenario 0 lives in slot 1.)

So it is a complete map of the dynamics — every arrow, precomputed.

Phase 2 then walks those arrows to find where each scenario ends up, without ever recomputing a score.

The two functions below do the same job as calling [`succession_step`](@ref) on every scenario in turn, and produce byte-identical results.
They exist because doing it the obvious way is far slower:

* succession_step allocates fresh vectors on each call (the scenario, its impact balance, the successor). Over hundreds of millions of scenarios that allocation dominates everything else.

* It also recomputes the impact balance from scratch each time, summing one matrix row per descriptor. But consecutive scenarios differ in only ONE descriptor, so almost all of that sum is unchanged from the previous scenario.

The version below therefore walks the scenarios in signature order like an odometer (see signatures.jl) and carries the impact balance along with it, changing it by the difference between two matrix rows whenever a digit ticks over. No allocation at all inside the loop.

"""

"""
    _successor_table!(successorTable, cib, cimTranspose) -> successorTable

Fill `successorTable[sig + 1]` with the successor of every scenario, as its signature.

The work is split into contiguous chunks of the scenario space, one task each. Because those ranges do not overlap, the tasks write to disjoint slots of the shared array and need no locking. Each chunk walks its own range with the odometer described above, and picks each descriptor's best variant using exactly the tie-breaking of [`succession_step`](@ref) under [`GlobalSuccession`](@ref) — ties favour the current variant, then the lowest index.
"""
function _successor_table!(successorTable::Vector{SignatureInt}, cib::CIB,
                           cimTranspose::Matrix{ScoreInt}) where {SignatureInt, ScoreInt}
    numberOfScenarios = max_signature(cib) + 1
    # The place value of each descriptor's digit in a signature — 1 for the
    # first descriptor, then multiplied by each variant count in turn (the
    # same mixed-radix weights `signature` uses). Precomputing them lets a
    # chunk build a successor's signature directly from the chosen variants,
    # so it never has to construct the successor scenario as a vector first.
    placeValues = Vector{Int}(undef, cib.numberOfDescriptors)
    placeValue = 1
    for descriptorIndex in 1:cib.numberOfDescriptors
        placeValues[descriptorIndex] = placeValue
        placeValue *= cib.numberOfVariants[descriptorIndex]
    end

    numberOfChunks = max(1, min(numberOfScenarios, 16 * Threads.nthreads()))
    chunkSize = cld(numberOfScenarios, numberOfChunks)
    numberOfChunks = cld(numberOfScenarios, chunkSize)
    @sync for chunkIndex in 1:numberOfChunks
        firstSignature = (chunkIndex - 1) * chunkSize
        lastSignature = min(chunkIndex * chunkSize, numberOfScenarios) - 1
        Threads.@spawn _successor_chunk!(successorTable, cimTranspose,
                                         firstSignature, lastSignature,
                                         cib.numberOfVariants, cib.desc_offsets,
                                         placeValues, cib.numberOfDescriptors,
                                         cib.numberOfDimensions)
    end
    return successorTable
end

# One worker's share of the table: fill in every scenario from firstSignature
# to lastSignature.
#
# This is the twin of _sweep_chunk_all! in sweep.jl — same odometer, same
# incremental impact balance. The only difference is what it does at each
# scenario. The sweep asks a yes/no question ("is anything better?") and can
# stop at the first descriptor that says yes; here we need every descriptor's
# actual choice, so there is no early exit.
function _successor_chunk!(successorTable::Vector{SignatureInt}, cimTranspose::Matrix{ScoreInt},
                           firstSignature::Int, lastSignature::Int,
                           variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                           placeValues::Vector{Int}, numberOfDescriptors::Int,
                           numberOfDimensions::Int) where {SignatureInt, ScoreInt}
    scenario      = Vector{Int}(undef, numberOfDescriptors)  # the odometer's current digits
    activeRows    = Vector{Int}(undef, numberOfDescriptors)  # matrix row of each descriptor's chosen variant
    impactBalance = zeros(ScoreInt, numberOfDimensions)      # carried along and patched, never rebuilt

    # Set the odometer to the chunk's starting scenario by splitting its
    # signature back into digits, and add up its impact balance from scratch.
    # This is the only division and the only full rescoring in the whole
    # chunk; from here on everything is incremental.
    remainder = firstSignature
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        variantCount = variantCounts[descriptorIndex]
        scenario[descriptorIndex] = remainder % variantCount
        activeRows[descriptorIndex] = descriptorOffsets[descriptorIndex] + scenario[descriptorIndex] + 1
        remainder = remainder ÷ variantCount
    end
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        sourceRow = activeRows[descriptorIndex]
        @simd for targetVariant in 1:numberOfDimensions
            impactBalance[targetVariant] += cimTranspose[targetVariant, sourceRow]
        end
    end

    @inbounds for currentSignature in firstSignature:lastSignature
        # ── Work out where this scenario steps to ──
        # Each descriptor picks its own best-scoring variant, exactly as the
        # global rule does. Rather than collecting those choices into a
        # scenario vector and converting it, we accumulate the signature
        # directly: each chosen variant contributes (its place value × its
        # number), which is precisely what `signature` computes.
        successorSignature = 0
        for descriptorIndex in 1:numberOfDescriptors
            offset = descriptorOffsets[descriptorIndex]
            # The variant currently in play is the incumbent, so ties keep it.
            bestVariant = scenario[descriptorIndex]
            bestScore = impactBalance[offset + bestVariant + 1]
            for variantIndex in 0:variantCounts[descriptorIndex]-1
                score = impactBalance[offset + variantIndex + 1]
                isBetter = score > bestScore   # strict >, as in succession_step
                # `ifelse` is not an `if`: it evaluates both alternatives and
                # then selects between them. That looks wasteful, but it
                # compiles to a conditional move rather than a jump, so the
                # CPU has no branch to mispredict. Worth it here because this
                # is the innermost loop of the whole analysis, and which
                # variant wins is essentially unpredictable.
                bestScore   = ifelse(isBetter, score, bestScore)
                bestVariant = ifelse(isBetter, variantIndex, bestVariant)
            end
            successorSignature += placeValues[descriptorIndex] * bestVariant
        end
        successorTable[currentSignature + 1] = SignatureInt(successorSignature)

        # ── Step the odometer on to the next scenario ──
        # Exactly the increment used in _sweep_chunk_all!: adding one changes
        # a single descriptor's variant, so the impact balance moves by the
        # difference between two matrix rows rather than being rebuilt.
        if currentSignature < lastSignature
            for descriptorIndex in 1:numberOfDescriptors
                variantCount = variantCounts[descriptorIndex]
                variantCount == 1 && continue    # single-variant digit never changes; carry onward
                oldRow = activeRows[descriptorIndex]
                if scenario[descriptorIndex] + 1 < variantCount
                    scenario[descriptorIndex] += 1
                    newRow = oldRow + 1
                    activeRows[descriptorIndex] = newRow
                    @simd for targetVariant in 1:numberOfDimensions
                        impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                        cimTranspose[targetVariant, oldRow]
                    end
                    break
                end
                scenario[descriptorIndex] = 0    # roll over; carry to the next digit
                newRow = descriptorOffsets[descriptorIndex] + 1
                activeRows[descriptorIndex] = newRow
                @simd for targetVariant in 1:numberOfDimensions
                    impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                    cimTranspose[targetVariant, oldRow]
                end
            end
        end
    end
    return successorTable
end

"""
    _fp_id!(registryLock, signatureForId, idForSignature, sig) -> id

Locked get-or-assign of a dense 1-based id for the fixed point with signature `sig`. Called once per fixed point discovered (≈ number-of-fixed-points times in total across all workers), so the lock is essentially uncontended. Concurrent discoverers of the same fixed point serialise here and receive the same id.
"""
# @noinline keeps this rarely-taken locked path out of the caller's hot
# loop, so the compiler optimises the loop without it.
@noinline function _fp_id!(registryLock::ReentrantLock, signatureForId::Vector{Int},
                           idForSignature::Dict{Int,Int}, fixedPointSignature::Int)
    # lock / try / finally-unlock is the standard exception-safe locking
    # pattern — the same shape as C#'s `lock` statement expands to.
    lock(registryLock)
    try
        denseId = get(idForSignature, fixedPointSignature, 0)  # 0 = not registered yet
        if denseId == 0
            # Labels are stored as Int32 (they hold dense ids, never signatures),
            # so more than ~2.1 billion DISTINCT fixed points cannot be labelled.
            # No real matrix is near this — it takes something degenerate like an
            # all-zero CIM over a >2^31 space, where every scenario is its own
            # fixed point — but wrapping silently would corrupt the tally, so:
            length(signatureForId) >= typemax(Int32) - 1 && throw(ArgumentError(
                "find_basins: more than $(typemax(Int32) - 1) distinct fixed points — " *
                "this degenerate model cannot be tallied with Int32 labels"))
            push!(signatureForId, fixedPointSignature)
            denseId = length(signatureForId)
            idForSignature[fixedPointSignature] = denseId
        end
        return denseId
    finally
        unlock(registryLock)
    end
end
#endregion

#region "Successors to consistent destination scenarios"
"""
    _resolve_chunk!(attractorLabels, successorTable, firstSignature, lastSignature,
                    registryLock, signatureForId, idForSignature)

Resolve the starting scenarios in `firstSignature:lastSignature` into the shared `attractorLabels`, following successor chains.
Cycle detection is **thread-private** (a per-worker `history` list plus a backward scan), so `attractorLabels` only ever holds *final* labels — 0 = unvisited, -1 = cycle, k > 0 = converges to the fixed point with dense id `k`. There is no shared in-progress marker, so a worker that walks into another worker's not-yet-resolved chain just re-walks it (redundant, never a false cycle) and reaches the same attractor; every scenario's attractor is deterministic, so concurrent writes to the same slot store the same value — a benign race. Fixed-point ids come from the locked registry ([`_fp_id!`](@ref)), hit only once per fixed point.

Labels are always `Int32`, whatever the width of the successor table's signatures: a label is a dense id or a sentinel, never a signature, and half-width labels halve the memory traffic of the resolve and tally passes (which are memory-bound pointer chases).
"""
function _resolve_chunk!(attractorLabels::Vector{Int32}, successorTable::Vector{SignatureInt},
                         firstSignature::Int, lastSignature::Int,
                         registryLock::ReentrantLock, signatureForId::Vector{Int},
                         idForSignature::Dict{Int,Int}) where {SignatureInt}
    history = Int[]     # the chain of scenarios walked from the current start
    @inbounds for startSignature in firstSignature:lastSignature
        attractorLabels[startSignature + 1] != 0 && continue   # already resolved
        empty!(history)
        current = startSignature
        label = zero(Int32)
        while true
            existingLabel = attractorLabels[current + 1]
            if existingLabel != 0
                # We walked into territory that is already resolved: the
                # whole chain behind us shares its attractor.
                label = existingLabel
                break
            end
            # Is `current` already on our own chain? Then we have walked in
            # a circle — a cycle of length ≥ 2. Scanning backwards finds a
            # repeat fastest, since a cycle closes near the chain's end.
            alreadyOnChain = false
            for position in length(history):-1:1
                if history[position] == current
                    alreadyOnChain = true
                    break
                end
            end
            if alreadyOnChain
                label = Int32(-1)   # the whole chain (incl. the tail leading in) is labelled cycle
                break
            end
            push!(history, current)
            next = Int(successorTable[current + 1])
            if next == current
                # A scenario that steps to itself is a fixed point. It is
                # already in `history`, so it counts itself in its own basin.
                label = Int32(_fp_id!(registryLock, signatureForId,
                                      idForSignature, current))
                break
            end
            current = next
        end
        # Path compression: every scenario on the walked chain gets the
        # final label, so later walks that touch any of them stop instantly.
        for visited in history
            attractorLabels[visited + 1] = label
        end
    end
    return nothing
end

"""
    _resolve_and_tally(successorTable, numberOfScenarios) -> (fp_sigs, sizes, cycle_count)

Resolve every scenario to its fixed point attractor by walking the successor table, then tally basin sizes and the cycle count.
"""
function _resolve_and_tally(successorTable::Vector{SignatureInt},
                            numberOfScenarios::Int) where {SignatureInt}
    # Always Int32: labels hold dense fixed-point ids (or 0/-1 sentinels), never
    # signatures, so they need no more width than the id space — see _fp_id!.
    attractorLabels = zeros(Int32, numberOfScenarios)
    registryLock = ReentrantLock()
    signatureForId = Int[]             # dense id (1-based) -> fixed-point signature
    idForSignature = Dict{Int,Int}()   # fixed-point signature -> dense id

    numberOfThreads = Threads.nthreads()
    chunkSize = cld(numberOfScenarios, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * chunkSize
        lastSignature = min(threadIndex * chunkSize, numberOfScenarios) - 1
        Threads.@spawn _resolve_chunk!(attractorLabels, successorTable,
                                       firstSignature, lastSignature,
                                       registryLock, signatureForId, idForSignature)
    end

    # ── Threaded tally. Every label is now a dense fixed-point id (> 0) or
    #    -1 for a cycle. Each thread counts into its own arrays; the counts
    #    are combined afterwards, so no locking is needed.
    numberOfFixedPoints = length(signatureForId)
    tallyChunkSize = cld(numberOfScenarios, numberOfThreads)
    perThreadCounts = [zeros(Int, numberOfFixedPoints) for _ in 1:numberOfThreads]
    perThreadCycleCounts = zeros(Int, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * tallyChunkSize
        lastSignature = min(threadIndex * tallyChunkSize, numberOfScenarios) - 1
        counts = perThreadCounts[threadIndex]
        Threads.@spawn begin
            cyclesSeen = 0
            @inbounds for signatureValue in firstSignature:lastSignature
                label = Int(attractorLabels[signatureValue + 1])
                if label == -1
                    cyclesSeen += 1
                else
                    counts[label] += 1   # label is a dense id in 1:numberOfFixedPoints
                end
            end
            perThreadCycleCounts[threadIndex] = cyclesSeen
        end
    end

    # Combine the per-thread counts into one total per fixed point.
    totalCounts = zeros(Int, numberOfFixedPoints)
    for counts in perThreadCounts
        @inbounds for fixedPointId in 1:numberOfFixedPoints
            totalCounts[fixedPointId] += counts[fixedPointId]
        end
    end
    # Discovery order depends on thread timing; sorting by signature makes
    # the output deterministic. `sortperm` returns the ordering as an index
    # list (like NumPy's argsort), applied to both arrays in step.
    sortOrder = sortperm(signatureForId)
    return signatureForId[sortOrder], totalCounts[sortOrder], sum(perThreadCycleCounts)
end
#endregion

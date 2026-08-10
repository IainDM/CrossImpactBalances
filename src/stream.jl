# ══ Streaming exact basins: every start, no tables ═════════════════════════
#
# find_basins' table method stores two flat arrays over the whole scenario
# space. At hundreds of millions of scenarios that is gigabytes and fine; at
# trillions it is tens of terabytes and impossible. This file computes the
# IDENTICAL result — exact basin sizes, exact cycle count — with memory that
# does not grow with the space at all: each worker walks every starting
# scenario in its share of the space to that trajectory's end (walk.jl) and
# keeps nothing but per-attractor counters.
#
# The price is honest arithmetic, not memory: with N starts and trajectories
# averaging L steps of ~ndesc·ndim work each, the stream does N·L step
# computations where the table did N. That is CPU-days at 10^12–10^13
# scenarios (split it across machines — scripts/basin_stream_worker.jl — and
# measure first with bench/stream_calibration.jl), and permanently out of
# reach at 10^18, where estimate_basins and influence_structure/product_basins
# are the tools.
#
# Three design points:
#
#   WORK SHARING. Trajectory lengths vary wildly between regions of the
#   space, and a run can last hours, so the fixed up-front chunking used
#   elsewhere in the codebase would strand threads behind one slow chunk.
#   Instead, workers CLAIM fixed-size blocks of starts off a shared atomic
#   counter as they finish — self-balancing, and still deterministic, because
#   per-attractor counts merge commutatively and the output is sorted by
#   fixed-point signature exactly as _resolve_and_tally does.
#
#   FIRST STEP FOR FREE. Starts are enumerated in signature order with the
#   same odometer as _successor_chunk!, carrying the start's impact balance
#   incrementally. So each start's FIRST successor costs one row-delta and an
#   argmax — only from the second step onwards does the walk pay full price.
#   Consistent scenarios (successor == start) are detected right there and
#   never walked at all.
#
#   THE MEMO CACHE. Trajectories funnel: different starts quickly merge onto
#   shared paths into the same attractor. A fixed-size, lossy, direct-mapped
#   cache of "signature → attractor id" lets a walk stop the moment it
#   touches territory any worker has already resolved. Each entry is ONE
#   UInt64 — the signature (as a tag) and the id packed together — so a
#   concurrent overwrite can never produce a half-entry made of two facts:
#   whatever whole word a reader sees is a true statement, and a mismatched
#   tag is simply a miss. (This extends the codebase's documented
#   benign-race pattern — see _resolve_chunk! — from "racing writers store
#   the same value" to "racing writers store different, self-validating
#   values". Aligned 8-byte stores do not tear on the platforms Julia
#   supports.) Collisions and evictions lose only speed, never correctness:
#   the walk that misses just keeps walking. The cache is used when the
#   packing fits (spaces up to 2^43 ≈ 8.8e12 scenarios) and skipped
#   otherwise; `cache_bytes=0` turns it off.

const _CACHE_ID_BITS = 21                              # up to ~2.1M distinct fixed-point ids
const _CACHE_ID_MASK = (UInt64(1) << _CACHE_ID_BITS) - 1
const _CACHE_CYCLE_ID = 0                              # id 0 = "this state reaches a cycle"
const _STREAM_BLOCK = 1 << 22                          # starts per claimed block (~4.2M)
const _STREAM_RING = 64                                # walk states remembered for insertion

struct _StreamCache
    slots::Vector{UInt64}
    mask::UInt64
end

"""
    _stream_cache(numberOfScenarios, cacheBytes) -> _StreamCache or nothing

Build the memo cache, or return `nothing` when it cannot help: when
`cacheBytes` is too small to matter, or when the space is so large that a
signature tag plus an id no longer pack into one 64-bit word (beyond 2^43
scenarios) — the stream then runs cacheless, exactly and more slowly.
"""
function _stream_cache(numberOfScenarios::Int, cacheBytes::Int)
    cacheBytes >= 8 * 1024 || return nothing
    # The tag is signature+1 (so an all-zero word means "empty slot"); it must
    # fit in the 64 - _CACHE_ID_BITS bits above the id.
    UInt64(numberOfScenarios) <= (UInt64(1) << (64 - _CACHE_ID_BITS)) - 1 || return nothing
    numberOfSlots = prevpow(2, max(1024, cacheBytes ÷ 8))
    return _StreamCache(zeros(UInt64, numberOfSlots), UInt64(numberOfSlots - 1))
end

# SplitMix64's finaliser as the slot hash: signatures are sequential, and a
# power-of-two table indexed by raw value would put a whole block of starts in
# adjacent slots; the mix scatters them.
@inline function _cache_slot(cache::_StreamCache, sig::Int)
    z = UInt64(sig)
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    z = z ⊻ (z >> 31)
    return (z & cache.mask) + 1
end

# -1 = miss; otherwise the stored attractor id (0 = cycle).
@inline function _cache_lookup(cache::_StreamCache, sig::Int)
    word = @inbounds cache.slots[_cache_slot(cache, sig)]
    word == 0 && return -1
    (word >> _CACHE_ID_BITS) == UInt64(sig) + 1 || return -1
    return Int(word & _CACHE_ID_MASK)
end

@inline function _cache_insert!(cache::_StreamCache, sig::Int, id::Int)
    # An id too large to pack is simply not cached (needs a pathological
    # million-fixed-point model); correctness never depends on an insert.
    id <= Int(_CACHE_ID_MASK) || return nothing
    @inbounds cache.slots[_cache_slot(cache, sig)] =
        ((UInt64(sig) + 1) << _CACHE_ID_BITS) | UInt64(id)
    return nothing
end
# No-op twins so the hot loop can call unconditionally when the cache is off.
@inline _cache_lookup(::Nothing, sig::Int) = -1
@inline _cache_insert!(::Nothing, sig::Int, id::Int) = nothing

"""
    _stream_walk!(space, ring, ringCount, cache, placeValues, ...) -> id

Walk from the scenario in `space.current` (a start's first successor) to
resolution, with three exits: a cache hit (some worker already resolved this
territory), a fixed point (registered through the shared [`_fp_id!`](@ref)
registry), or a Brent-detected cycle. Every state visited has the same fate
as the start, so their signatures are remembered in the fixed-size `ring`
(most recent `_STREAM_RING` of them — the tail end of a walk is the shared
funnel, which is exactly the part worth caching) and inserted into the cache
with the resolved id before returning.

Returns the attractor id (0 = cycle) — the caller tallies it.
"""
function _stream_walk!(space::_WalkSpace{ScoreInt}, ring::Vector{Int}, ringCount::Int,
                       cache::Union{_StreamCache,Nothing},
                       placeValues::Vector{Int}, cimTranspose::Matrix{ScoreInt},
                       variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                       numberOfDescriptors::Int, numberOfDimensions::Int,
                       registryLock::ReentrantLock, signatureForId::Vector{Int},
                       idForSignature::Dict{Int,Int}) where {ScoreInt<:Signed}
    copyto!(space.anchor, space.current)     # Brent's saved scenario (see walk.jl)
    power = 1
    stretch = 0
    id = _CACHE_CYCLE_ID
    while true
        becameFixed = _global_step!(space.scratch, space.current, space.balance,
                                    cimTranspose, variantCounts, descriptorOffsets,
                                    numberOfDescriptors, numberOfDimensions)
        space.current, space.scratch = space.scratch, space.current
        stretch += 1

        # The new state's signature, from the digits (no per-step allocation).
        sig = 0
        @inbounds for descriptorIndex in 1:numberOfDescriptors
            sig += placeValues[descriptorIndex] * space.current[descriptorIndex]
        end

        if becameFixed
            id = _fp_id!(registryLock, signatureForId, idForSignature, sig)
            ringCount = _ring_push!(ring, ringCount, sig)
            break
        end
        cached = _cache_lookup(cache, sig)
        if cached >= 0
            id = cached
            break                     # sig itself is already cached; no need to re-add
        end
        ringCount = _ring_push!(ring, ringCount, sig)

        if space.current == space.anchor
            id = _CACHE_CYCLE_ID      # walked in a circle: a cycle of length `stretch`
            break
        end
        if power == stretch
            copyto!(space.anchor, space.current)
            power <<= 1
            stretch = 0
        end
    end

    # Everything this walk visited shares the resolved fate: publish it.
    for position in 1:min(ringCount, _STREAM_RING)
        _cache_insert!(cache, ring[position], id)
    end
    return id
end

@inline function _ring_push!(ring::Vector{Int}, ringCount::Int, sig::Int)
    @inbounds ring[(ringCount % _STREAM_RING) + 1] = sig
    return ringCount + 1
end

"""
    _stream_block!(...) -> (cyclesSeen, stepsWalked)

One claimed block of starts: enumerate `firstSignature:lastSignature` with the
incremental odometer (the same preamble and tick as `_successor_chunk!`),
resolve each start via first-step / cache / walk, and tally into the worker's
private `tally` (attractor id → count). Returns this block's cycle count and
total succession steps (the latter feeds calibration statistics).
"""
function _stream_block!(tally::Dict{Int,Int}, space::_WalkSpace{ScoreInt}, ring::Vector{Int},
                        cache::Union{_StreamCache,Nothing}, cimTranspose::Matrix{ScoreInt},
                        firstSignature::Int, lastSignature::Int,
                        variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                        placeValues::Vector{Int}, numberOfDescriptors::Int,
                        numberOfDimensions::Int, registryLock::ReentrantLock,
                        signatureForId::Vector{Int},
                        idForSignature::Dict{Int,Int}) where {ScoreInt<:Signed}
    scenario      = Vector{Int}(undef, numberOfDescriptors)
    activeRows    = Vector{Int}(undef, numberOfDescriptors)
    impactBalance = zeros(ScoreInt, numberOfDimensions)
    cyclesSeen = 0
    stepsWalked = 0

    # Odometer preamble: decode the block's first start, score it from scratch.
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
        # ── The start's first successor, from the carried balance ──
        # The same incumbent-seeded strict-> argmax as everywhere else, but
        # accumulating both the successor's signature and its digits (into
        # space.current, where a walk would begin).
        successorSignature = 0
        anythingMoved = false
        for descriptorIndex in 1:numberOfDescriptors
            offset = descriptorOffsets[descriptorIndex]
            bestVariant = scenario[descriptorIndex]
            bestScore = impactBalance[offset + bestVariant + 1]
            for variantIndex in 0:variantCounts[descriptorIndex]-1
                score = impactBalance[offset + variantIndex + 1]
                isBetter = score > bestScore
                bestScore   = ifelse(isBetter, score, bestScore)
                bestVariant = ifelse(isBetter, variantIndex, bestVariant)
            end
            space.current[descriptorIndex] = bestVariant
            successorSignature += placeValues[descriptorIndex] * bestVariant
            anythingMoved |= (bestVariant != scenario[descriptorIndex])
        end
        stepsWalked += 1

        if !anythingMoved
            # The start is itself a consistent scenario: its own basin member.
            id = _fp_id!(registryLock, signatureForId, idForSignature, currentSignature)
            tally[id] = get(tally, id, 0) + 1
            _cache_insert!(cache, currentSignature, id)
        else
            id = _cache_lookup(cache, successorSignature)
            if id < 0
                # Prime the ring with the two states already known to share
                # the walk's fate, then walk from the successor.
                ringCount = _ring_push!(ring, 0, currentSignature)
                ringCount = _ring_push!(ring, ringCount, successorSignature)
                id = _stream_walk!(space, ring, ringCount, cache, placeValues,
                                   cimTranspose, variantCounts, descriptorOffsets,
                                   numberOfDescriptors, numberOfDimensions,
                                   registryLock, signatureForId, idForSignature)
            else
                # Successor already resolved: the start inherits its fate.
                _cache_insert!(cache, currentSignature, id)
            end
            if id == _CACHE_CYCLE_ID
                cyclesSeen += 1
            else
                tally[id] = get(tally, id, 0) + 1
            end
        end

        # ── Odometer tick (identical to _successor_chunk!) ──
        if currentSignature < lastSignature
            for descriptorIndex in 1:numberOfDescriptors
                variantCount = variantCounts[descriptorIndex]
                variantCount == 1 && continue
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
                scenario[descriptorIndex] = 0
                newRow = descriptorOffsets[descriptorIndex] + 1
                activeRows[descriptorIndex] = newRow
                @simd for targetVariant in 1:numberOfDimensions
                    impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                    cimTranspose[targetVariant, oldRow]
                end
            end
        end
    end
    return cyclesSeen, stepsWalked
end

# One worker: claim blocks off the shared counter until none remain.
function _stream_worker!(tally::Dict{Int,Int}, cib::CIB, cache::Union{_StreamCache,Nothing},
                         cimTranspose::Matrix{ScoreInt}, nextBlock::Threads.Atomic{Int},
                         blocksDone::Threads.Atomic{Int}, numberOfBlocks::Int,
                         firstSignature::Int, lastSignature::Int, placeValues::Vector{Int},
                         registryLock::ReentrantLock, signatureForId::Vector{Int},
                         idForSignature::Dict{Int,Int}) where {ScoreInt<:Signed}
    space = _WalkSpace(cib, ScoreInt)
    ring = Vector{Int}(undef, _STREAM_RING)
    cyclesSeen = 0
    stepsWalked = 0
    while true
        blockIndex = Threads.atomic_add!(nextBlock, 1)   # returns the value BEFORE the add
        blockIndex >= numberOfBlocks && break
        blockFirst = firstSignature + blockIndex * _STREAM_BLOCK
        blockLast = min(blockFirst + _STREAM_BLOCK - 1, lastSignature)
        cycles, steps = _stream_block!(tally, space, ring, cache, cimTranspose,
                                       blockFirst, blockLast, cib.numberOfVariants,
                                       cib.desc_offsets, placeValues,
                                       cib.numberOfDescriptors, cib.numberOfDimensions,
                                       registryLock, signatureForId, idForSignature)
        cyclesSeen += cycles
        stepsWalked += steps
        Threads.atomic_add!(blocksDone, 1)
    end
    return cyclesSeen, stepsWalked
end

"""
    _stream_basins(rule, cib, firstSignature, lastSignature, cacheBytes, progress)
        -> (fixedPoints, basinSizes, cycleCount)

The streaming counterpart of the table method, over the starts in
`firstSignature:lastSignature` — the whole space unless [`find_basins`](@ref)
was given a `signature_range` to split the work across machines. Output is
identical in form and (over the full space) in content to the table method's,
including the ascending-fixed-point-signature ordering and the
`sum(sizes) + cycles == number of starts` invariant.
"""
function _stream_basins(::GlobalSuccession, cib::CIB, firstSignature::Int, lastSignature::Int,
                        cacheBytes::Int, progress::Bool)
    scoreType = _score_type(cib)
    return _stream_basins_fast(cib, Matrix{scoreType}(cib.cim_t),
                               firstSignature, lastSignature, cacheBytes, progress)
end

function _stream_basins_fast(cib::CIB, cimTranspose::Matrix{ScoreInt},
                             firstSignature::Int, lastSignature::Int,
                             cacheBytes::Int, progress::Bool) where {ScoreInt<:Signed}
    numberOfScenarios = _basin_scenario_count(cib)
    placeValues = Vector{Int}(undef, cib.numberOfDescriptors)
    placeValue = 1
    for descriptorIndex in 1:cib.numberOfDescriptors
        placeValues[descriptorIndex] = placeValue
        placeValue *= cib.numberOfVariants[descriptorIndex]
    end

    cache = cacheBytes > 0 ? _stream_cache(numberOfScenarios, cacheBytes) : nothing
    registryLock = ReentrantLock()
    signatureForId = Int[]
    idForSignature = Dict{Int,Int}()

    numberOfStarts = lastSignature - firstSignature + 1
    numberOfBlocks = cld(numberOfStarts, _STREAM_BLOCK)
    numberOfWorkers = max(1, min(Threads.nthreads(), numberOfBlocks))
    perWorkerTallies = [Dict{Int,Int}() for _ in 1:numberOfWorkers]
    perWorkerCycles = zeros(Int, numberOfWorkers)
    perWorkerSteps = zeros(Int, numberOfWorkers)
    nextBlock = Threads.Atomic{Int}(0)
    blocksDone = Threads.Atomic{Int}(0)

    monitor = progress ? _stream_monitor(blocksDone, numberOfBlocks, numberOfStarts) : nothing

    @sync for workerIndex in 1:numberOfWorkers
        Threads.@spawn begin
            cycles, steps = _stream_worker!(perWorkerTallies[workerIndex], cib, cache,
                                            cimTranspose, nextBlock, blocksDone,
                                            numberOfBlocks, firstSignature, lastSignature,
                                            placeValues, registryLock, signatureForId,
                                            idForSignature)
            perWorkerCycles[workerIndex] = cycles
            perWorkerSteps[workerIndex] = steps
        end
    end
    monitor === nothing || close(monitor)

    # Merge the private tallies. Integer addition commutes, so the merge order
    # (and therefore the thread count and block schedule) cannot affect totals.
    totalCounts = zeros(Int, length(signatureForId))
    for tally in perWorkerTallies
        for (denseId, count) in tally
            totalCounts[denseId] += count
        end
    end
    if progress
        walked = sum(perWorkerSteps)
        println(stderr, "find_basins stream: done — $numberOfStarts starts, ",
                walked, " succession steps (", round(walked / numberOfStarts, digits=2),
                " per start incl. cache effects)")
    end

    sortOrder = sortperm(signatureForId)               # deterministic output order
    fixedPoints = [inv_signature(cib, sig) for sig in signatureForId[sortOrder]]
    return (fixedPoints, totalCounts[sortOrder], sum(perWorkerCycles))
end

# A Timer prints progress every few seconds from the scheduler, without a
# dedicated thread; close()ing it stops the printing. Plain stderr lines, the
# same convention the MCP worker and app use.
function _stream_monitor(blocksDone::Threads.Atomic{Int}, numberOfBlocks::Int,
                         numberOfStarts::Int)
    startedAt = time()
    return Timer(5.0; interval=5.0) do _
        done = blocksDone[]
        elapsed = time() - startedAt
        rate = done * _STREAM_BLOCK / max(elapsed, 1e-9)
        remaining = rate > 0 ? (numberOfStarts - done * _STREAM_BLOCK) / rate : NaN
        println(stderr, "find_basins stream: ", done, "/", numberOfBlocks,
                " blocks  (~", round(Int, rate), " starts/s, ~",
                round(remaining / 3600, digits=1), " h remaining)")
    end
end

# ── Generic-rule stream ─────────────────────────────────────────────────────
# Any custom SuccessionRule gets the same streaming analysis through its own
# succession_step — the exact mirror of find_basins' fast/generic split. No
# memo cache here: the rule's step allocates anyway, and correctness needs
# nothing from the cache.
function _stream_basins(rule::SuccessionRule, cib::CIB, firstSignature::Int, lastSignature::Int,
                        cacheBytes::Int, progress::Bool)
    registryLock = ReentrantLock()
    signatureForId = Int[]
    idForSignature = Dict{Int,Int}()

    numberOfStarts = lastSignature - firstSignature + 1
    numberOfBlocks = cld(numberOfStarts, _STREAM_BLOCK)
    numberOfWorkers = max(1, min(Threads.nthreads(), numberOfBlocks))
    perWorkerTallies = [Dict{Int,Int}() for _ in 1:numberOfWorkers]
    perWorkerCycles = zeros(Int, numberOfWorkers)
    nextBlock = Threads.Atomic{Int}(0)
    blocksDone = Threads.Atomic{Int}(0)
    monitor = progress ? _stream_monitor(blocksDone, numberOfBlocks, numberOfStarts) : nothing

    @sync for workerIndex in 1:numberOfWorkers
        Threads.@spawn begin
            tally = perWorkerTallies[workerIndex]
            cyclesSeen = 0
            while true
                blockIndex = Threads.atomic_add!(nextBlock, 1)
                blockIndex >= numberOfBlocks && break
                blockFirst = firstSignature + blockIndex * _STREAM_BLOCK
                blockLast = min(blockFirst + _STREAM_BLOCK - 1, lastSignature)
                for startSignature in blockFirst:blockLast
                    start = inv_signature(cib, startSignature)
                    isFixedPoint, attractor, _ = _walk_to_attractor(rule, cib, start)
                    if isFixedPoint
                        denseId = _fp_id!(registryLock, signatureForId, idForSignature,
                                          signature(cib, attractor))
                        tally[denseId] = get(tally, denseId, 0) + 1
                    else
                        cyclesSeen += 1
                    end
                end
                Threads.atomic_add!(blocksDone, 1)
            end
            perWorkerCycles[workerIndex] = cyclesSeen
        end
    end
    monitor === nothing || close(monitor)

    totalCounts = zeros(Int, length(signatureForId))
    for tally in perWorkerTallies
        for (denseId, count) in tally
            totalCounts[denseId] += count
        end
    end
    sortOrder = sortperm(signatureForId)
    fixedPoints = [inv_signature(cib, sig) for sig in signatureForId[sortOrder]]
    return (fixedPoints, totalCounts[sortOrder], sum(perWorkerCycles))
end

# ══ Walking a single succession trajectory, without a table ════════════════
#
# find_basins' fast path precomputes every scenario's successor into one huge
# flat table and then chases pointers through it. That is the right design
# when the table fits in memory — each succession is computed exactly once —
# and impossible when it doesn't: at 6.7 trillion scenarios the table alone
# would need tens of terabytes.
#
# This file is the other way to follow the dynamics: start somewhere, compute
# one successor at a time, and keep going until the trajectory ends. No table,
# no memory proportional to the scenario count — just three scenario-sized
# buffers and a score vector. Both the sampling analysis (estimate.jl) and the
# streaming exact analysis (stream.jl) are built on it.
#
# Two things need care:
#
#   ENDING. A trajectory in a finite space must eventually repeat — either a
#   fixed point (a consistent scenario, repeating with period 1) or a longer
#   cycle. We detect the repeat with BRENT'S ALGORITHM: keep one saved
#   "anchor" scenario, and compare the moving scenario against it; every time
#   the walk doubles in length, move the anchor up to the current position.
#   The anchor is guaranteed to eventually sit ON the cycle while the walker
#   comes around again, so a repeat is always caught — using two buffers and
#   zero growing state, unlike the visited-set approach elsewhere in the
#   codebase, which this analysis cannot afford billions of times over.
#
#   SPEED. Each step recomputes the impact balance in full: ndesc row-sums of
#   the transposed matrix (the SIMD-friendly layout impact_balance uses), then
#   each descriptor picks its best variant under exactly the tie rule of
#   succession_step — the incumbent opens the contest and only a STRICTLY
#   better challenger unseats it. That inner argmax is deliberately written
#   out again here rather than shared with _successor_chunk! (basins.jl) or
#   the sweep: those are tuned hot loops the codebase keeps as twins on
#   purpose, and equivalence is pinned by the property tests instead.

"""
    _WalkSpace(cib, ScoreInt)

The preallocated working memory for walking trajectories: three scenario
buffers (`current`, `anchor`, `scratch`) and one score accumulator. One
_WalkSpace per worker task; every walk reuses it, so a walk allocates nothing.
`mutable` so the buffers can be swapped by reference instead of copied.
"""
mutable struct _WalkSpace{ScoreInt<:Signed}
    current::Vector{Int}    # the walker: where the trajectory is now
    anchor::Vector{Int}     # Brent's saved scenario, compared against current
    scratch::Vector{Int}    # the successor under construction / cycle minimum
    balance::Vector{ScoreInt}
end

_WalkSpace(cib::CIB, ::Type{ScoreInt}) where {ScoreInt<:Signed} =
    _WalkSpace{ScoreInt}(Vector{Int}(undef, cib.numberOfDescriptors),
                         Vector{Int}(undef, cib.numberOfDescriptors),
                         Vector{Int}(undef, cib.numberOfDescriptors),
                         zeros(ScoreInt, cib.numberOfDimensions))

"""
    _global_step!(destination, source, balance, cimTranspose, variantCounts,
                  descriptorOffsets, numberOfDescriptors, numberOfDimensions) -> Bool

One step of [`GlobalSuccession`](@ref) computed from scratch: fill `balance`
with the impact balance of `source`, then write each descriptor's winning
variant into `destination`. Returns `true` when nothing moved — i.e. `source`
is a consistent scenario. Produces byte-identical choices to
[`succession_step`](@ref) (the property tests hold the two together).
"""
@inline function _global_step!(destination::Vector{Int}, source::Vector{Int},
                               balance::Vector{ScoreInt}, cimTranspose::Matrix{ScoreInt},
                               variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                               numberOfDescriptors::Int, numberOfDimensions::Int) where {ScoreInt<:Signed}
    fill!(balance, zero(ScoreInt))
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        sourceRow = descriptorOffsets[descriptorIndex] + source[descriptorIndex] + 1
        @simd for targetVariant in 1:numberOfDimensions
            balance[targetVariant] += cimTranspose[targetVariant, sourceRow]
        end
    end

    anythingMoved = false
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        offset = descriptorOffsets[descriptorIndex]
        # The variant currently in play opens the contest, so ties keep it —
        # the same strict-> rule as succession_step and the odometer paths.
        bestVariant = source[descriptorIndex]
        bestScore = balance[offset + bestVariant + 1]
        for variantIndex in 0:variantCounts[descriptorIndex]-1
            score = balance[offset + variantIndex + 1]
            isBetter = score > bestScore
            bestScore   = ifelse(isBetter, score, bestScore)
            bestVariant = ifelse(isBetter, variantIndex, bestVariant)
        end
        destination[descriptorIndex] = bestVariant
        anythingMoved |= (bestVariant != source[descriptorIndex])
    end
    return !anythingMoved
end

"""
    _sig_less(a, b) -> Bool

Is scenario `a`'s signature smaller than scenario `b`'s? Decided without any
signature arithmetic: signatures are mixed-radix numbers whose MOST
significant digit is the LAST descriptor, so comparing the digits from the
last descriptor downwards is exactly numeric signature order — and it works
even for models whose signatures would overflow `Int64`.
"""
@inline function _sig_less(a::Vector{Int}, b::Vector{Int})
    @inbounds for descriptorIndex in length(a):-1:1
        a[descriptorIndex] != b[descriptorIndex] &&
            return a[descriptorIndex] < b[descriptorIndex]
    end
    return false
end

"""
    _walk_to_attractor!(space, start, cimTranspose, variantCounts,
                        descriptorOffsets, numberOfDescriptors, numberOfDimensions)
        -> (isFixedPoint, steps)

Follow [`GlobalSuccession`](@ref) from `start` until the trajectory's end,
using Brent's cycle detection (see the file header). On return the attractor
is in `space.current`: the consistent scenario reached (`isFixedPoint ==
true`), or — for a cycle — the cycle's smallest-signature member, a canonical
choice so that every walk entering the same cycle names it identically.
`steps` counts succession steps taken, for calibration statistics. `start` is
not modified.
"""
function _walk_to_attractor!(space::_WalkSpace{ScoreInt}, start::Vector{Int},
                             cimTranspose::Matrix{ScoreInt},
                             variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                             numberOfDescriptors::Int, numberOfDimensions::Int) where {ScoreInt<:Signed}
    copyto!(space.anchor, start)
    # current = successor(start). If nothing moved, start was already consistent.
    if _global_step!(space.current, space.anchor, space.balance, cimTranspose,
                     variantCounts, descriptorOffsets, numberOfDescriptors, numberOfDimensions)
        return (true, 1)
    end
    steps = 1

    # ── Brent's search for the first repeat ──
    # `power` is the doubling schedule; `stretch` counts steps since the anchor
    # was last moved. When the walker equals the anchor, `stretch` is the cycle
    # length. A fixed point shows up as a cycle of length 1.
    power = 1
    stretch = 1
    while space.current != space.anchor
        if power == stretch
            copyto!(space.anchor, space.current)
            power <<= 1
            stretch = 0
        end
        # scratch = successor(current), then swap the two buffers by reference.
        _global_step!(space.scratch, space.current, space.balance, cimTranspose,
                      variantCounts, descriptorOffsets, numberOfDescriptors, numberOfDimensions)
        space.current, space.scratch = space.scratch, space.current
        steps += 1
        stretch += 1
    end

    cycleLength = stretch
    cycleLength == 1 && return (true, steps)   # period 1: a consistent scenario

    # ── A genuine cycle: pick its canonical (smallest-signature) member ──
    # Walk once around the loop from where we stand, keeping the minimum seen.
    copyto!(space.anchor, space.current)       # the running minimum
    for _ in 1:cycleLength-1
        _global_step!(space.scratch, space.current, space.balance, cimTranspose,
                      variantCounts, descriptorOffsets, numberOfDescriptors, numberOfDimensions)
        space.current, space.scratch = space.scratch, space.current
        steps += 1
        _sig_less(space.current, space.anchor) && copyto!(space.anchor, space.current)
    end
    copyto!(space.current, space.anchor)       # contract: attractor ends up in `current`
    return (false, steps)
end

"""
    _walk_to_attractor(rule, cib, start) -> (isFixedPoint, attractor, steps)

The generic-rule twin of [`_walk_to_attractor!`](@ref): the same Brent walk,
driven by the rule's own [`succession_step`](@ref) so it works for ANY
[`SuccessionRule`](@ref) — at the cost of the allocations succession_step
makes per step. Mirrors the fast-path/generic-path split of
[`find_basins`](@ref) itself.
"""
function _walk_to_attractor(rule::SuccessionRule, cib::CIB, start::Vector{Int})
    anchor = copy(start)
    current = succession_step(rule, cib, anchor)
    current == anchor && return (true, current, 1)
    steps = 1

    power = 1
    stretch = 1
    while current != anchor
        if power == stretch
            anchor = current
            power <<= 1
            stretch = 0
        end
        current = succession_step(rule, cib, current)
        steps += 1
        stretch += 1
    end

    cycleLength = stretch
    cycleLength == 1 && return (true, current, steps)

    minimum = current
    for _ in 1:cycleLength-1
        current = succession_step(rule, cib, current)
        steps += 1
        _sig_less(current, minimum) && (minimum = current)
    end
    return (false, minimum, steps)
end

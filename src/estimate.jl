# ══ Estimating basin shares by sampling, with honest error bars ═════════════
#
# Beyond ~10^13 scenarios, no exact method can visit every start: the table
# won't fit and even the memory-flat stream would need CPU-centuries. What CAN
# be had at any scale, in seconds, is an estimate of each basin's SHARE of the
# space, with confidence intervals: draw starting scenarios uniformly at
# random, walk each to its end (walk.jl), and count where they land. A basin
# holding 31% of a 10^18 space announces itself after a few thousand draws;
# the error bars shrink with 1/√samples and do not care how big the space is.
#
# THE DIVISION OF LABOUR MATTERS, and it is what keeps this package's promise
# intact. Finding the consistent scenarios stays EXACT — find_consistent's
# branch-and-bound already handles these spaces, and every consistent scenario
# is registered in the result before the first sample is drawn. Sampling is
# never asked to DISCOVER anything, only to apportion mass among attractors
# already known — so a tiny basin cannot silently vanish the way it does when
# sampling does the discovering (the documented CIBSA failure on bench_50x50:
# a 10k-point sample missing 2 of 5 consistent scenarios, whose basins were 3
# and 95 out of 60.5M). Here a never-hit attractor still appears, with an
# explicit upper bound on its share instead of a false zero.
#
# Determinism: the estimate is a pure function of (model, rule, samples, seed,
# confidence) — NOT of the thread count. Samples are partitioned into
# fixed-size blocks by sample index alone; each block draws from its own
# generator derived from (seed, block index) alone (rng.jl); tallies are
# integer counts merged commutatively and sorted before assembly. Same seed,
# same answer, any machine, any -t.

const _ESTIMATE_BLOCK = 1 << 16    # samples per block: the unit of reproducibility

"""
    BasinEstimate

The result of [`estimate_basins`](@ref): estimated basin shares with
confidence intervals, for a kernel that is known exactly.

Parallel vectors, sorted by ascending attractor signature:
- `fixedPoints`: every attractor in the report — the exact kernel supplied or
  computed up front, plus any fixed point the sampler reached (there should be
  none of the latter if the kernel was complete)
- `hits`: how many sampled starts drained into each (0 is meaningful: the
  attractor exists, and its share was too small to see at this sample size)
- `shares`: `hits / samples`
- `ciLow`, `ciHigh`: the Wilson interval for each share at level `confidence`
- `sizeEstimates`: `shares × scenarioCount`, as `Float64` — at these scales a
  size estimate is meaningful to its leading digits only

And the remainder of the tally:
- `cycleHits`, `cycleShare`, `cycleCiLow`, `cycleCiHigh`: starts that ended in
  a repeating cycle (no basin), as count and interval
- `samples`, `seed`, `confidence`: the estimate's defining inputs
- `scenarioCount::Int128`: the exact size of the space ([`scenario_count`](@ref))

`sum(hits) + cycleHits == samples` always holds — the sampled twin of the
exact analysis' `sum(sizes) + cycles == total` invariant.
"""
struct BasinEstimate
    fixedPoints::Vector{Vector{Int}}
    hits::Vector{Int}
    shares::Vector{Float64}
    ciLow::Vector{Float64}
    ciHigh::Vector{Float64}
    sizeEstimates::Vector{Float64}
    cycleHits::Int
    cycleShare::Float64
    cycleCiLow::Float64
    cycleCiHigh::Float64
    samples::Int
    seed::UInt64
    confidence::Float64
    scenarioCount::Int128
end

"""
    estimate_basins(cib; rule=GlobalSuccession(), samples=1_000_000,
                    seed=0xC1BBA512C1BBA512, confidence=0.95, kernel=nothing)
        -> BasinEstimate

Estimate each consistent scenario's basin SHARE by uniform random sampling —
the analysis for scenario spaces beyond every exact method's reach (they run
out around 10¹³; this does not care). Each of `samples` starting scenarios is
drawn uniformly (every scenario in the space equally likely), walked to its
end under `rule`, and tallied; shares come back with Wilson confidence
intervals at level `confidence`.

**This estimates; it never discovers.** The kernel is taken exact — from
`kernel` if given, else the model's already-computed consistent scenarios,
else [`find_consistent`](@ref) (whose branch-and-bound handles huge spaces) —
and every consistent scenario appears in the result even with zero hits,
carrying an explicit upper bound on its share rather than a false zero. A
share the sampler cannot see (much below ~1/`samples`) is reported exactly
that way. For rules that declare no [`fixed_point_margin`](@ref) and supply no
`kernel`, the kernel cannot be precomputed at this scale; the estimate then
reports only what sampling reaches, and says so.

The result is a pure function of `(model, rule, samples, seed, confidence)` —
independent of thread count, machine, and Julia version (the generator is the
package's own; see rng.jl). Same seed, same estimate, bit for bit.

See [`find_basins`](@ref) for exact analysis when the space allows it, and
[`influence_structure`](@ref) / [`product_basins`](@ref) for exact analysis of
huge models that decompose.
"""
function estimate_basins(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                         samples::Integer=1_000_000,
                         seed::UInt64=0xC1BBA512C1BBA512,
                         confidence::Float64=0.95,
                         kernel::Union{Nothing,Vector{Vector{Int}}}=nothing)
    samples >= 1 || throw(ArgumentError("estimate_basins: samples must be ≥ 1, got $samples"))
    z = _z_for_confidence(confidence)          # validates confidence ∈ (0, 1)
    totalScenarios = scenario_count(cib)

    kernelScenarios, kernelKnown = _resolve_kernel(cib, rule, kernel)

    # ── Shared tally state ──
    sampleCount = Int(samples)
    numberOfBlocks = cld(sampleCount, _ESTIMATE_BLOCK)
    numberOfWorkers = max(1, min(Threads.nthreads(), numberOfBlocks))
    perWorkerTallies = [Dict{Int128,Int}() for _ in 1:numberOfWorkers]
    perWorkerScenarios = [Dict{Int128,Vector{Int}}() for _ in 1:numberOfWorkers]
    perWorkerCycles = zeros(Int, numberOfWorkers)
    nextBlock = Threads.Atomic{Int}(0)

    scoreType = _score_type(cib)
    cimTranspose = rule isa GlobalSuccession ? Matrix{scoreType}(cib.cim_t) :
                                               Matrix{scoreType}(undef, 0, 0)

    @sync for workerIndex in 1:numberOfWorkers
        Threads.@spawn begin
            perWorkerCycles[workerIndex] =
                _estimate_worker!(perWorkerTallies[workerIndex],
                                  perWorkerScenarios[workerIndex], rule, cib,
                                  cimTranspose, nextBlock, numberOfBlocks,
                                  sampleCount, seed)
        end
    end

    # ── Merge (commutative: integer sums; first-writer scenario copies are
    #    identical by construction, since a signature determines its scenario) ──
    mergedHits = Dict{Int128,Int}()
    mergedScenarios = Dict{Int128,Vector{Int}}()
    for workerIndex in 1:numberOfWorkers
        for (key, count) in perWorkerTallies[workerIndex]
            mergedHits[key] = get(mergedHits, key, 0) + count
        end
        for (key, scenario) in perWorkerScenarios[workerIndex]
            haskey(mergedScenarios, key) || (mergedScenarios[key] = scenario)
        end
    end
    # Pre-registered kernel: present in the report even with zero hits.
    for scenario in kernelScenarios
        key = _signature128(cib, scenario)
        haskey(mergedHits, key) || (mergedHits[key] = 0)
        haskey(mergedScenarios, key) || (mergedScenarios[key] = copy(scenario))
    end

    sortedKeys = sort!(collect(keys(mergedHits)))
    hits = [mergedHits[key] for key in sortedKeys]
    fixedPoints = [mergedScenarios[key] for key in sortedKeys]
    shares = [h / sampleCount for h in hits]
    intervals = [_wilson(h, sampleCount, z) for h in hits]
    cycleHits = sum(perWorkerCycles)
    cycleInterval = _wilson(cycleHits, sampleCount, z)

    return BasinEstimate(fixedPoints, hits, shares,
                         [interval[1] for interval in intervals],
                         [interval[2] for interval in intervals],
                         [share * Float64(totalScenarios) for share in shares],
                         cycleHits, cycleHits / sampleCount,
                         cycleInterval[1], cycleInterval[2],
                         sampleCount, seed, confidence, totalScenarios)
end

# Decide where the exact kernel comes from. Returns (scenarios, known::Bool);
# known == false means no kernel could be established (margin-less custom rule,
# none supplied) and the report can only show what sampling happened to reach.
function _resolve_kernel(cib::CIB, rule::SuccessionRule,
                         kernel::Union{Nothing,Vector{Vector{Int}}})
    resolved = if kernel !== nothing
        kernel
    elseif fixed_point_margin(rule) == 0 && !isempty(cib.consistentScenarios)
        # Any margin-0 rule has exactly the model's consistent scenarios as its
        # fixed points (that is what declaring margin 0 promises), so the
        # already-computed kernel is reusable as-is.
        cib.consistentScenarios
    elseif fixed_point_margin(rule) !== nothing
        find_consistent(cib; rule=rule)
    else
        println(stderr, "estimate_basins: this rule declares no fixed_point_margin and no " *
                        "kernel was supplied — the exact kernel cannot be precomputed at " *
                        "this scale, so the report covers only attractors that sampling " *
                        "reaches. Pass kernel=... to pre-register known fixed points.")
        return Vector{Vector{Int}}(), false
    end
    for scenario in resolved
        length(scenario) == cib.numberOfDescriptors || throw(ArgumentError(
            "estimate_basins: kernel scenario $scenario has $(length(scenario)) descriptors, " *
            "model has $(cib.numberOfDescriptors)"))
        succession_step(rule, cib, scenario) == scenario || throw(ArgumentError(
            "estimate_basins: kernel scenario $scenario is not a fixed point under this rule"))
    end
    return resolved, true
end

# One worker: claim sample blocks until none remain. The block index alone
# determines the generator and which sample indices the block covers, so the
# assignment of blocks to workers is irrelevant to the result.
function _estimate_worker!(tally::Dict{Int128,Int}, scenarios::Dict{Int128,Vector{Int}},
                           rule::SuccessionRule, cib::CIB,
                           cimTranspose::Matrix{ScoreInt}, nextBlock::Threads.Atomic{Int},
                           numberOfBlocks::Int, sampleCount::Int,
                           seed::UInt64) where {ScoreInt<:Signed}
    space = _WalkSpace(cib, ScoreInt)
    start = Vector{Int}(undef, cib.numberOfDescriptors)
    cyclesSeen = 0
    while true
        blockIndex = Threads.atomic_add!(nextBlock, 1)
        blockIndex >= numberOfBlocks && break
        rng = _block_rng(seed, blockIndex)
        blockSamples = min(_ESTIMATE_BLOCK, sampleCount - blockIndex * _ESTIMATE_BLOCK)
        cyclesSeen += _estimate_block!(tally, scenarios, rule, cib, cimTranspose,
                                       space, start, rng, blockSamples)
    end
    return cyclesSeen
end

# Fast path: GlobalSuccession through the preallocated walker.
function _estimate_block!(tally::Dict{Int128,Int}, scenarios::Dict{Int128,Vector{Int}},
                          ::GlobalSuccession, cib::CIB, cimTranspose::Matrix{ScoreInt},
                          space::_WalkSpace{ScoreInt}, start::Vector{Int},
                          rng::_Xoshiro256pp, blockSamples::Int) where {ScoreInt<:Signed}
    cyclesSeen = 0
    for _ in 1:blockSamples
        @inbounds for descriptorIndex in 1:cib.numberOfDescriptors
            start[descriptorIndex] =
                Int(_rand_below!(rng, UInt64(cib.numberOfVariants[descriptorIndex])))
        end
        isFixedPoint, _ = _walk_to_attractor!(space, start, cimTranspose,
                                              cib.numberOfVariants, cib.desc_offsets,
                                              cib.numberOfDescriptors, cib.numberOfDimensions)
        if isFixedPoint
            key = _signature128(cib, space.current)
            tally[key] = get(tally, key, 0) + 1
            haskey(scenarios, key) || (scenarios[key] = copy(space.current))
        else
            cyclesSeen += 1
        end
    end
    return cyclesSeen
end

# Generic path: any rule, through its own succession_step.
function _estimate_block!(tally::Dict{Int128,Int}, scenarios::Dict{Int128,Vector{Int}},
                          rule::SuccessionRule, cib::CIB, ::Matrix{ScoreInt},
                          ::_WalkSpace{ScoreInt}, start::Vector{Int},
                          rng::_Xoshiro256pp, blockSamples::Int) where {ScoreInt<:Signed}
    cyclesSeen = 0
    for _ in 1:blockSamples
        @inbounds for descriptorIndex in 1:cib.numberOfDescriptors
            start[descriptorIndex] =
                Int(_rand_below!(rng, UInt64(cib.numberOfVariants[descriptorIndex])))
        end
        isFixedPoint, attractor, _ = _walk_to_attractor(rule, cib, start)
        if isFixedPoint
            key = _signature128(cib, attractor)
            tally[key] = get(tally, key, 0) + 1
            haskey(scenarios, key) || (scenarios[key] = attractor)
        else
            cyclesSeen += 1
        end
    end
    return cyclesSeen
end

# ── Display ─────────────────────────────────────────────────────────────────

# 1234567 -> "1,234,567" (also for Int128); purely cosmetic.
function _group_digits(n::Integer)
    digits = string(abs(widen(n)))
    grouped = replace(reverse(join(Iterators.partition(reverse(digits), 3), ",")), r"^," => "")
    return (n < 0 ? "-" : "") * grouped
end

_pct(x::Float64) = string(round(100 * x, sigdigits=4)) * "%"

function Base.show(io::IO, ::MIME"text/plain", estimate::BasinEstimate)
    println(io, "BasinEstimate: ", _group_digits(estimate.samples), " samples over ",
            _group_digits(estimate.scenarioCount), " scenarios")
    println(io, "  seed 0x", string(estimate.seed, base=16), ", ",
            _pct(estimate.confidence), " Wilson intervals")
    # Largest share first for reading; the struct's vectors stay signature-sorted.
    order = sortperm(estimate.hits; rev=true)
    for index in order
        scenario = estimate.fixedPoints[index]
        if estimate.hits[index] == 0
            println(io, "  ", scenario, ": no hits — share ≤ ",
                    _pct(estimate.ciHigh[index]), " (upper bound)")
        else
            println(io, "  ", scenario, ": ", _pct(estimate.shares[index]),
                    "  [", _pct(estimate.ciLow[index]), ", ", _pct(estimate.ciHigh[index]),
                    "]  ≈ ", round(estimate.sizeEstimates[index], sigdigits=3), " scenarios")
        end
    end
    println(io, "  cycles (no basin): ", _pct(estimate.cycleShare),
            "  [", _pct(estimate.cycleCiLow), ", ", _pct(estimate.cycleCiHigh), "]")
    print(io, "  Exact kernel; estimated shares — discovery is never delegated to sampling.")
end

Base.show(io::IO, estimate::BasinEstimate) =
    print(io, "BasinEstimate(", length(estimate.fixedPoints), " attractors, ",
          _group_digits(estimate.samples), " samples)")

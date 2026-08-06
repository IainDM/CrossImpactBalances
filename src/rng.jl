# ══ Random numbers and interval statistics, self-contained ═════════════════
#
# The sampling analysis (estimate.jl) needs a stream of random numbers that is
# REPRODUCIBLE: the same seed must give bit-identical estimates on any machine,
# any Julia version, any thread count, forever — the same determinism contract
# the exact analyses honour (bench/determinism.jl). Julia's own `Random` module
# deliberately does NOT promise stream stability across minor releases, and
# this package ships with no dependencies at all, stdlib included. So the ~40
# lines of generator are written out here instead.
#
# The algorithms are the standard, public-domain ones (Blackman & Vigna):
#
#   SPLITMIX64     A tiny generator whose only job here is turning one 64-bit
#                  seed into several well-scrambled state words.
#   XOSHIRO256++   The main generator — the same family Julia's own default
#                  RNG uses. Fast, tiny state, passes the standard test
#                  batteries; entirely adequate for Monte Carlo tallying.
#
# Plus two pieces of interval statistics for the estimator's error bars: an
# inverse normal CDF (Acklam's rational approximation) and the Wilson score
# interval for a binomial proportion.

"""
    _splitmix64(state) -> (output, newstate)

One step of the SplitMix64 generator. Used only for seeding: it turns
consecutive states into decorrelated 64-bit outputs, which become the initial
state words of [`_Xoshiro256pp`](@ref). (Returning a tuple is Julia's normal
way to hand back two values; there are no out-parameters.)
"""
@inline function _splitmix64(state::UInt64)
    state += 0x9E3779B97F4A7C15          # golden-ratio increment; wraps mod 2^64 by design
    z = state
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9   # ⊻ is xor
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    return z ⊻ (z >> 31), state
end

"""
    _Xoshiro256pp(seed::UInt64)

The xoshiro256++ generator, seeded by expanding `seed` through four rounds of
[`_splitmix64`](@ref) (the seeding procedure its authors recommend). A
`mutable struct` because each call to [`_rand64!`](@ref) advances the state
in place.
"""
mutable struct _Xoshiro256pp
    s0::UInt64
    s1::UInt64
    s2::UInt64
    s3::UInt64
end

function _Xoshiro256pp(seed::UInt64)
    word0, state = _splitmix64(seed)
    word1, state = _splitmix64(state)
    word2, state = _splitmix64(state)
    word3, _     = _splitmix64(state)
    # The all-zero state is xoshiro's one forbidden state (it would emit zeros
    # forever). SplitMix64 producing four zero words is not going to happen by
    # chance, but the guard costs nothing and removes the caveat.
    if word0 == word1 == word2 == word3 == 0
        word0 = 0x9E3779B97F4A7C15
    end
    return _Xoshiro256pp(word0, word1, word2, word3)
end

# Rotate the bits of x left by k — single instruction on every modern CPU.
@inline _rotl(x::UInt64, k::Int) = (x << k) | (x >> (64 - k))

"""
    _rand64!(rng) -> UInt64

The next 64 uniformly random bits, advancing the generator state.
"""
@inline function _rand64!(rng::_Xoshiro256pp)
    result = _rotl(rng.s0 + rng.s3, 23) + rng.s0
    t = rng.s1 << 17
    rng.s2 ⊻= rng.s0
    rng.s3 ⊻= rng.s1
    rng.s1 ⊻= rng.s2
    rng.s0 ⊻= rng.s3
    rng.s2 ⊻= t
    rng.s3 = _rotl(rng.s3, 45)
    return result
end

"""
    _block_rng(seed, blockIndex) -> _Xoshiro256pp

An independent generator for one block of samples, derived from the master
`seed` and the block's index — and from nothing else. This is what makes the
sampling analysis independent of the thread count: block 17 draws the same
numbers whether it runs on thread 1 of 1 or thread 8 of 8, so the estimate is
a pure function of (seed, sample count). The golden-ratio multiply spreads
consecutive block indices far apart in seed space before the SplitMix64
expansion scrambles them.
"""
_block_rng(seed::UInt64, blockIndex::Integer) =
    _Xoshiro256pp(seed ⊻ (UInt64(blockIndex) * 0x9E3779B97F4A7C15))

"""
    _rand_below!(rng, bound) -> UInt64

A uniform integer in `0:bound-1`, unbiased for every bound ≥ 1 (Lemire's
multiply-shift method, with the rejection step). The naive `x % bound` is
imperceptibly biased; an estimator whose whole point is calibrated error bars
should not be built on a biased die, however slightly. The occasional
rejection re-draw is part of the deterministic stream — same seed, same
rejections, same result. For CIB variant counts (2–9ish) a rejection almost
never happens.
"""
@inline function _rand_below!(rng::_Xoshiro256pp, bound::UInt64)
    x = _rand64!(rng)
    m = UInt128(x) * UInt128(bound)      # widen: the product needs 128 bits
    low = m % UInt64                      # low 64 bits of the product
    if low < bound
        threshold = (-bound) % bound      # = 2^64 mod bound, computed in UInt64
        while low < threshold
            x = _rand64!(rng)
            m = UInt128(x) * UInt128(bound)
            low = m % UInt64
        end
    end
    return (m >> 64) % UInt64             # high 64 bits: the unbiased result
end

# ── Interval statistics ─────────────────────────────────────────────────────

"""
    _inv_normal_cdf(p) -> Float64

The standard normal quantile Φ⁻¹(p) by Acklam's rational approximation
(relative error < 1.2e-9 — far below anything a Monte Carlo error bar can
resolve). Self-contained so the estimator does not need a statistics package.
"""
function _inv_normal_cdf(p::Float64)
    0.0 < p < 1.0 || throw(ArgumentError("_inv_normal_cdf: p must be in (0, 1), got $p"))
    # Break points between the central rational fit and the two tail fits.
    pLow = 0.02425
    if p < pLow                                        # lower tail
        q = sqrt(-2.0 * log(p))
        return (((((-7.784894002430293e-3 * q - 3.223964580411365e-1) * q -
                   2.400758277161838e0) * q - 2.549732539343734e0) * q +
                 4.374664141464968e0) * q + 2.938163982698783e0) /
               ((((7.784695709041462e-3 * q + 3.224671290700398e-1) * q +
                  2.445134137142996e0) * q + 3.754408661907416e0) * q + 1.0)
    elseif p <= 1.0 - pLow                             # central region
        q = p - 0.5
        r = q * q
        return (((((-3.969683028665376e1 * r + 2.209460984245205e2) * r -
                   2.759285104469687e2) * r + 1.383577518672690e2) * r -
                 3.066479806614716e1) * r + 2.506628277459239e0) * q /
               (((((-5.447609879822406e1 * r + 1.615858368580409e2) * r -
                   1.556989798598866e2) * r + 6.680131188771972e1) * r -
                 1.328068155288572e1) * r + 1.0)
    else                                               # upper tail (mirror of lower)
        q = sqrt(-2.0 * log(1.0 - p))
        return -(((((-7.784894002430293e-3 * q - 3.223964580411365e-1) * q -
                    2.400758277161838e0) * q - 2.549732539343734e0) * q +
                  4.374664141464968e0) * q + 2.938163982698783e0) /
                ((((7.784695709041462e-3 * q + 3.224671290700398e-1) * q +
                   2.445134137142996e0) * q + 3.754408661907416e0) * q + 1.0)
    end
end

"""
    _z_for_confidence(confidence) -> Float64

The two-sided normal critical value for a confidence level, e.g. ≈1.96 for
0.95.
"""
function _z_for_confidence(confidence::Float64)
    0.0 < confidence < 1.0 ||
        throw(ArgumentError("confidence must be in (0, 1), got $confidence"))
    return _inv_normal_cdf(1.0 - (1.0 - confidence) / 2.0)
end

"""
    _wilson(hits, n, z) -> (low, high)

The Wilson score interval for a binomial proportion with `hits` successes in
`n` trials at critical value `z`. Chosen over the naive ±z·√(p(1-p)/n)
interval because it behaves sensibly at the edges this analysis actually
lives at: `hits == 0` gives a positive upper bound rather than a zero-width
interval (that upper bound, ≈ z²/n ≈ 3.84/n at 95%, is how a basin the
sampler never hit gets an honest "no larger than..." statement — the
rule-of-three idea, done properly), and small counts don't produce intervals
poking below 0 or above 1.
"""
function _wilson(hits::Int, n::Int, z::Float64)
    n > 0 || throw(ArgumentError("_wilson: need n > 0, got $n"))
    0 <= hits <= n || throw(ArgumentError("_wilson: hits=$hits outside 0:$n"))
    observedShare = hits / n
    zSquaredOverN = z^2 / n
    denominator = 1.0 + zSquaredOverN
    center = (observedShare + zSquaredOverN / 2.0) / denominator
    halfWidth = z * sqrt(observedShare * (1.0 - observedShare) / n +
                         zSquaredOverN / (4.0 * n)) / denominator
    return (max(0.0, center - halfWidth), min(1.0, center + halfWidth))
end

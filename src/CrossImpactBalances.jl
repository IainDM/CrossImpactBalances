"""
    CrossImpactBalances

Cross-Impact Balance (CIB) scenario analysis. A Julia implementation of the
methodology in Weimer-Jehle (2006), ported from the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with a two-phase basin-of-attraction analysis, a threaded incremental
exhaustive sweep, and an exact branch-and-bound search that prunes
provably-inconsistent regions of the scenario space. The succession dynamics
is pluggable via [`SuccessionRule`](@ref): [`GlobalSuccession`](@ref) (the
default, carrying the fast paths) and any user-defined rule drop into the
same analysis routines through multiple dispatch.

# Reference
Weimer-Jehle, W. (2006). Cross-impact balances: A system-theoretical approach
to cross-impact analysis. *Technological Forecasting and Social Change*,
73(4), 334–361.
"""
module CrossImpactBalances

using SIMD: Vec, vload, vifelse

export CIB, load_scw, load_solutions,
       impact_balance,
       SuccessionRule, GlobalSuccession, SequentialSuccession,
       succession_step, succession,
       find_consistent, find_basins,
       signature, inv_signature, max_signature

"""
    CIB

Cross-impact balance analysis object.

Fields:
- `descriptors`: ordered list of descriptor names
- `variants`: dict mapping descriptor name to list of variant names
- `nvariants`: number of variants per descriptor
- `cim`: the cross-impact matrix (n × n), n = sum of all variants. Row `r`
  holds the impacts contributed by variant `r` onto every other variant.
- `cim_t`: the transpose of `cim` — `cim_t[j, r] == cim[r, j]`. Stored
  alongside `cim` so that "sum a fixed-row vector across many `j`" hot
  loops (impact_balance, find_basins) can iterate contiguously through
  column-major storage, which unlocks SIMD. Construction cost is one
  `permutedims` at load time.
- `ndim`: total number of variants (size of CIM)
- `ndesc`: number of descriptors
- `kernel`: list of consistent scenarios (each a Vector{Int}, 0-based variant indices)
- `desc_offsets`: cumulative 0-based variant offsets per descriptor

`CIB` is immutable in its field bindings. `kernel` is a `Vector` whose
contents may be mutated in place.
"""
struct CIB
    descriptors::Vector{String}
    variants::Dict{String, Vector{String}}
    nvariants::Vector{Int}
    cim::Matrix{Int}
    cim_t::Matrix{Int}
    ndim::Int
    ndesc::Int
    kernel::Vector{Vector{Int}}
    desc_offsets::Vector{Int}
end

# ─── .scw file parser ───────────────────────────────────────────────────────

"""
    load_scw(scw_file; sl_file=nothing, kernel=nothing, algorithm=:auto) -> CIB

Parse a ScenarioWizard .scw file and optionally a .sl solutions file.
Returns a fully populated CIB object.

Unless a `kernel` or an `sl_file` is supplied, the kernel is computed by
[`find_consistent`](@ref), which always searches the full scenario space and
uses every available thread (start Julia with `julia -t auto`). `algorithm`
is forwarded to select the search strategy.
"""
function load_scw(scw_file::String; sl_file::Union{String,Nothing}=nothing,
                  kernel::Union{Vector{Vector{Int}},Nothing}=nothing,
                  algorithm::Symbol=:auto)
    descriptors = String[]
    variants = Dict{String, Vector{String}}()
    nvars = Int[]

    # State machine
    s = 0
    n = 0       # total variant count
    d = -1      # descriptor counter
    v = 0       # variant counter for current descriptor
    r = 0       # row counter for CIM
    current_desc = ""

    # First pass: count dimensions
    cim_rows = Vector{Vector{Int}}()

    for line in eachline(scw_file)
        stripped = lstrip(line)
        isempty(stripped) && continue

        if s == 0
            if stripped[1] == '&'
                desc = strip(stripped[2:end])
                push!(descriptors, desc)
                variants[desc] = String[]
                if d > -1
                    push!(nvars, v)
                end
                v = 0
                d += 1
                current_desc = desc
            elseif stripped[1] == '-'
                push!(variants[current_desc], strip(stripped[2:end]))
                n += 1
                v += 1
            elseif stripped[1] == '#'
                s = 1
                push!(nvars, v)
            end
        elseif s < 5 && stripped[1] == '#'
            s += 1
        elseif s == 5
            if stripped[1] == '#'
                s += 1
            else
                row = parse.(Int, split(stripped, ','))
                push!(cim_rows, row)
                r += 1
            end
        end
    end

    n == 0 && error("load_scw: no variants found in $(scw_file) — file is empty or malformed")
    length(cim_rows) == n || error("load_scw: cross-impact matrix has $(length(cim_rows)) " *
                                   "rows but $n variants in $(scw_file)")

    # Build CIM matrix
    cim = zeros(Int, n, n)
    for (i, row) in enumerate(cim_rows)
        length(row) == n || error("load_scw: row $i of CIM has $(length(row)) entries but $n expected")
        for (j, val) in enumerate(row)
            cim[i, j] = val
        end
    end

    ndesc = d + 1

    desc_offsets = Vector{Int}(undef, ndesc)
    off = 0
    for i in 1:ndesc
        desc_offsets[i] = off
        off += nvars[i]
    end

    # Precompute the transpose so impact_balance / find_basins can do
    # contiguous SIMD column reads in cim_t (= the row vectors of cim).
    cim_t = permutedims(cim)

    # Build CIB without a kernel first; if needed, populate kernel in place.
    cib = CIB(descriptors, variants, nvars, cim, cim_t, n, ndesc,
              Vector{Vector{Int}}(), desc_offsets)

    if !isnothing(kernel)
        append!(cib.kernel, kernel)
    elseif !isnothing(sl_file)
        append!(cib.kernel, load_solutions(cib, sl_file))
    else
        append!(cib.kernel, find_consistent(cib; algorithm=algorithm))
    end

    return cib
end

# ─── .sl file parser ────────────────────────────────────────────────────────

"""
    load_solutions(cib::CIB, sl_file::String) -> Vector{Vector{Int}}

Parse a ScenarioWizard .sl solutions file. Returns 0-based variant indices.
"""
function load_solutions(cib::CIB, sl_file::String)
    kern = Vector{Vector{Int}}()
    for line in eachline(sl_file)
        stripped = lstrip(line)
        isempty(stripped) && continue
        stripped[1] != '"' && continue

        # Extract the quoted index string, e.g. "2 3 2"
        m = match(r"^\"([^\"]+)\"", stripped)
        isnothing(m) && continue

        indices = parse.(Int, split(strip(m.captures[1])))
        # Convert from 1-based (ScenarioWizard) to 0-based (internal)
        push!(kern, indices .- 1)
    end
    return kern
end

# ─── Signature (unique integer ID for a scenario) ──────────────────────────

"""
    signature(cib, u) -> Int

Compute a unique integer signature for scenario `u` (0-based variant indices).
"""
function signature(cib::CIB, u::Vector{Int})
    sig = 0
    order = 1
    for (ui, nv) in zip(u, cib.nvariants)
        sig += order * ui
        order *= nv
    end
    return sig
end

"""
    inv_signature(cib, s) -> Vector{Int}

Convert a signature back to a scenario (0-based variant indices).
"""
function inv_signature(cib::CIB, s::Int)
    u = Int[]
    for nv in cib.nvariants
        push!(u, s % nv)
        s = s ÷ nv
    end
    return u
end

"""
    max_signature(cib) -> Int

The maximum signature value (= total scenarios - 1).
"""
function max_signature(cib::CIB)
    return signature(cib, [nv - 1 for nv in cib.nvariants])
end

# ─── Impact balance ─────────────────────────────────────────────────────────

"""
    impact_balance(cib, u) -> Vector{Int}

Compute the impact balance vector for scenario `u`.
Returns a vector of length ndim (one score per variant across all descriptors).
"""
function impact_balance(cib::CIB, u::Vector{Int})
    ib = zeros(Int, cib.ndim)
    nd = cib.ndim
    cim_t = cib.cim_t                      # row r of cim lives at column r of cim_t
    offsets = cib.desc_offsets
    @inbounds for (i, ui) in enumerate(u)
        r = offsets[i] + ui + 1
        @simd for j in 1:nd
            ib[j] += cim_t[j, r]           # contiguous column read → AVX-512
        end
    end
    return ib
end


# ─── Succession rules (pluggable dynamics) ──────────────────────────────────

"""
    SuccessionRule

Abstract supertype for a *succession dynamics* — the deterministic map that
sends a scenario to its successor. A consistent scenario is a fixed point of
this map; a basin of attraction is the set of scenarios that flow to a given
fixed point under repeated application.

Swapping in a new dynamics is the whole extension interface. Define

    struct MyRule <: SuccessionRule end

and a single method

    succession_step(rule::MyRule, cib::CIB, u::Vector{Int}) -> Vector{Int}

and every analysis routine — [`succession`](@ref), [`find_consistent`](@ref)
and [`find_basins`](@ref) — works with it immediately through a generic,
rule-agnostic path. [`find_basins`](@ref) is already fast for any rule (it
fills a successor table and reuses the shared path-compressed resolver). A
*threshold* rule — one whose fixed points are a separable, impact-balance-only
per-descriptor condition — can additionally opt [`find_consistent`](@ref) into
the fast threaded sweep and branch-and-bound by defining
[`fixed_point_margin`](@ref); otherwise `find_consistent` uses a generic scan.
Both are optional — correctness never depends on them.
"""
abstract type SuccessionRule end

"""
    GlobalSuccession()

Global (simultaneous) succession: in one step, every descriptor is moved to
the variant with the highest impact score given the *current* scenario. Ties
favour the current variant, then the lowest index. This is the classical
ScenarioWizard / CIBSA dynamics and the package default; it carries the fast
threaded sweep, branch-and-bound, and two-phase basin implementations.
"""
struct GlobalSuccession <: SuccessionRule end

"""
    SequentialSuccession()

Sequential (successive / Gauss–Seidel) succession: descriptors are updated one
at a time in descriptor order, each using the impact balance of the scenario
*as already partially updated within the same step*. An alternative CIB
dynamics, included mainly to demonstrate that the analysis routines are
rule-agnostic — it plugs into [`find_consistent`](@ref) and
[`find_basins`](@ref) through the generic path with no engine changes.
"""
struct SequentialSuccession <: SuccessionRule end

"""
    fixed_point_margin(rule) -> Union{Nothing,Int}

Fast-path opt-in for a *threshold* succession rule. Return an integer margin
`m ≥ 0` when this rule's fixed points are exactly the scenarios in which, for
every descriptor, no competing variant's impact score exceeds the current
variant's by more than `m` — i.e. the fixed-point test is separable per
descriptor and a function of the impact balance alone. Return `nothing` (the
default) for any rule without that structure.

Declaring a margin lets [`find_consistent`](@ref) reach the fast threaded
sweep and branch-and-bound (the margin parameterises their per-descriptor
fixed-point test) instead of the generic per-scenario scan. The two must
agree, so a margin is a *promise* about the rule's fixed points worth
validating in tests. [`GlobalSuccession`](@ref) is the `m = 0` case: a
competitor strictly beating the current variant unseats it.

Note this concerns *fixed points* only, not the full dynamics — a rule may be
a threshold rule for consistency yet have basins that still need its own
[`succession_step`](@ref) (e.g. sequential/Gauss–Seidel updates).
"""
fixed_point_margin(::SuccessionRule) = nothing
fixed_point_margin(::GlobalSuccession) = 0

# ─── Succession ─────────────────────────────────────────────────────────────

"""
    succession_step(cib, u) -> Vector{Int}
    succession_step(rule, cib, u) -> Vector{Int}

One step of succession under `rule` (default [`GlobalSuccession`](@ref)).
For global succession, each descriptor independently picks the variant with
the highest impact score given the current scenario; ties favour the current
variant, then the lowest index.
"""
succession_step(cib::CIB, u::Vector{Int}) = succession_step(GlobalSuccession(), cib, u)

function succession_step(::GlobalSuccession, cib::CIB, u::Vector{Int})
    ib = impact_balance(cib, u)
    v = copy(u)
    start = 1  # 1-based index into ib
    for i in 1:cib.ndesc
        nv = cib.nvariants[i]
        stop = start + nv - 1
        ib_desc = @view ib[start:stop]
        max_val = ib_desc[u[i] + 1]  # current variant's score (+1 for 1-based)
        for j in 0:nv-1
            if ib_desc[j + 1] > max_val  # strict >, so ties keep current/lower index
                max_val = ib_desc[j + 1]
                v[i] = j
            end
        end
        start = stop + 1
    end
    return v
end

function succession_step(::SequentialSuccession, cib::CIB, u::Vector{Int})
    v = copy(u)
    @inbounds for i in 1:cib.ndesc
        ib = impact_balance(cib, v)     # recomputed from the partially-updated v
        off = cib.desc_offsets[i]
        nv = cib.nvariants[i]
        max_val = ib[off + v[i] + 1]    # current variant's score (favour on ties)
        for j in 0:nv-1
            if ib[off + j + 1] > max_val
                max_val = ib[off + j + 1]
                v[i] = j
            end
        end
    end
    return v
end

"""
    succession(cib, u) -> (cycle_length, attractor)
    succession(rule, cib, u) -> (cycle_length, attractor)

Follow succession under `rule` (default [`GlobalSuccession`](@ref)) from
scenario `u` until convergence to a fixed point or detection of a cycle.
Returns (cycle_length, final_scenario). cycle_length=1 means a consistent
scenario (fixed point).

Cycle detection uses an O(1)-per-step hashtable, so total work is linear in
the trajectory length.
"""
succession(cib::CIB, u::Vector{Int}) = succession(GlobalSuccession(), cib, u)

function succession(rule::SuccessionRule, cib::CIB, u::Vector{Int})
    start_sig = signature(cib, u)
    history_sig = Int[start_sig]
    seen = Dict{Int,Int}()  # signature -> 1-based position in history_sig
    seen[start_sig] = 1
    v = copy(u)
    while true
        v = succession_step(rule, cib, v)
        v_sig = signature(cib, v)
        if haskey(seen, v_sig)
            cycle_len = length(history_sig) - seen[v_sig] + 1
            return (cycle_len, v)
        end
        push!(history_sig, v_sig)
        seen[v_sig] = length(history_sig)
    end
end

# ─── Find consistent scenarios ──────────────────────────────────────────────

"""
    find_consistent(cib; rule=GlobalSuccession(), algorithm=:auto,
                    bnb_node_budget=nothing) -> Vector{Vector{Int}}

Find every consistent scenario — every fixed point of the succession map under
`rule` — by exhaustively searching the whole scenario space with all available
threads (start Julia with `julia -t auto`).

`rule` selects the succession *dynamics* ([`SuccessionRule`](@ref); default
[`GlobalSuccession`](@ref)). A rule that declares a [`fixed_point_margin`](@ref)
— including the default `GlobalSuccession` (`m = 0`) and any threshold/inertial
rule — uses the fast sweep / branch-and-bound kernel; any other rule uses a
generic ascending-signature scan.

`algorithm` selects the search strategy (only for rules with a
[`fixed_point_margin`](@ref)):
- `:sweep` — enumerate every scenario with the incremental odometer sweep.
- `:bnb`   — branch-and-bound: assign descriptors depth-first and prune
  subtrees that provably contain no fixed point.
  Exact — returns the identical kernel — and typically visits a small
  fraction of the space on strongly-coupled matrices.
- `:auto`  (default) — the sweep for small spaces (< 10^5 scenarios),
  otherwise branch-and-bound with a node budget of `n ÷ 16`; if pruning is
  too weak to pay off, the budget trips and the sweep runs instead.

The returned kernel is ordered by ascending signature.
"""
function find_consistent(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                         algorithm::Symbol=:auto,
                         bnb_node_budget::Union{Nothing,Int}=nothing)
    return _exhaustive_kernel(rule, cib; algorithm=algorithm,
                              bnb_node_budget=bnb_node_budget)
end

"""
    _exhaustive_kernel(rule, cib; algorithm, bnb_node_budget)

Exhaustive fixed-point search under `rule`. A rule that declares a
[`fixed_point_margin`](@ref) gets the fast threaded sweep / branch-and-bound —
the margin parameterises the per-descriptor fixed-point test, so both paths
serve global succession (`m = 0`), inertial/threshold rules (`m > 0`) and any
other rule with the same separable structure. A rule with no margin falls back
to a generic ascending-signature scan that tests `succession_step(rule, u) == u`.
"""
function _exhaustive_kernel(rule::SuccessionRule, cib::CIB; algorithm::Symbol=:auto,
                            bnb_node_budget::Union{Nothing,Int}=nothing)
    margin = fixed_point_margin(rule)
    if margin === nothing
        # Generic fallback: enumerate the whole space in ascending-signature
        # order and keep the scenarios that are their own successor. Single-
        # threaded O(n) correctness baseline for a rule with no separable path.
        algorithm === :auto || throw(ArgumentError(
            "algorithm=$(repr(algorithm)) needs a rule with a fixed_point_margin; " *
            "this rule uses the generic scan (algorithm=:auto)"))
        kern = Vector{Vector{Int}}()
        for sig in 0:max_signature(cib)
            u = inv_signature(cib, sig)
            succession_step(rule, cib, u) == u && push!(kern, u)
        end
        return kern
    end
    # Threshold rule: the fast paths, parameterised by the margin.
    algorithm in (:auto, :bnb, :sweep) ||
        throw(ArgumentError("algorithm must be :auto, :bnb or :sweep, got $(repr(algorithm))"))
    n = max_signature(cib) + 1
    if algorithm == :sweep || (algorithm == :auto && n < 100_000)
        return _find_consistent_exhaustive(cib; margin=margin)
    end
    sufmin, sufmax = _bnb_bounds(cib)
    budget = something(bnb_node_budget,
                       algorithm == :bnb ? typemax(Int) : n ÷ 16)
    kern, _ = _bnb_fixed_points(cib, sufmin, sufmax; node_budget=budget, margin=margin)
    !isnothing(kern) && return kern
    # Budget tripped: pruning too weak on this matrix — fall back.
    return _find_consistent_exhaustive(cib; margin=margin)
end

"""
    _score_type(cib) -> Type{<:Signed}

Narrowest safe integer type for impact-balance accumulation. `Int16` when
`4 * (ndesc + 1) * maximum(abs, cim)` fits (margin covers the incremental
row-delta updates; Int16 arithmetic wraps silently, so this guard is the
only protection), otherwise `Int`. Realistic CIB matrices (entries ±3)
are far inside the Int16 bound.
"""
function _score_type(cib::CIB)
    maxabs = cib.ndim == 0 ? 0 : Int(maximum(abs, cib.cim))
    return 4 * (cib.ndesc + 1) * maxabs <= Int(typemax(Int16)) ? Int16 : Int
end

"""
    _find_consistent_exhaustive(cib) -> Vector{Vector{Int}}

Internal threaded exhaustive search. Splits the signature space into
`16 × Threads.nthreads()` contiguous chunks scheduled as tasks (early-exit
cost varies across the space, so fine-grained chunks load-balance better
than one block per thread). Within a chunk, the impact-balance vector is
maintained *incrementally*: a mixed-radix odometer advances the scenario,
and each digit change updates the balance by the difference of two CIM
rows (a `@simd` loop over the narrowed transpose from [`_score_type`](@ref)),
so no per-scenario score recomputation happens at all. The fixed-point
check is then a per-descriptor comparison over the balance vector with
early exit.

The returned kernel is ordered by ascending signature. Only true fixed points
(cycle length 1) are emitted, so non-fixed-point cycles never appear.
"""
function _find_consistent_exhaustive(cib::CIB; margin::Int=0)
    T = _score_type(cib)
    return _sweep_fixed_points(cib, Matrix{T}(cib.cim_t); margin=margin)
end

function _sweep_fixed_points(cib::CIB, cimT::Matrix{T}; margin::Int=0) where {T<:Signed}
    n = max_signature(cib) + 1
    nchunks = max(1, min(n, 16 * Threads.nthreads()))
    chunk_size = cld(n, nchunks)
    nchunks = cld(n, chunk_size)

    results = [Vector{Vector{Int}}() for _ in 1:nchunks]
    @sync for c in 1:nchunks
        out = results[c]
        first_sig = (c - 1) * chunk_size
        last_sig = min(c * chunk_size, n) - 1
        Threads.@spawn _sweep_chunk!(out, cimT, first_sig, last_sig,
                                     cib.nvariants, cib.desc_offsets,
                                     cib.ndesc, cib.ndim, margin)
    end

    # Chunks cover ascending contiguous signature ranges; merging in chunk
    # order keeps the kernel sorted by signature. No dedup needed.
    kern = Vector{Vector{Int}}()
    for c in 1:nchunks
        append!(kern, results[c])
    end
    return kern
end

function _sweep_chunk!(out::Vector{Vector{Int}}, cimT::Matrix{T},
                       first_sig::Int, last_sig::Int,
                       nvariants::Vector{Int}, offsets::Vector{Int},
                       ndesc::Int, ndim::Int, margin::Int) where {T<:Signed}
    v    = Vector{Int}(undef, ndesc)
    rows = Vector{Int}(undef, ndesc)   # column of cimT holding descriptor i's current row
    ib   = zeros(T, ndim)

    # Decode first_sig (the only divmod in this chunk) and build the
    # initial impact balance from scratch.
    s = first_sig
    @inbounds for i in 1:ndesc
        nv = nvariants[i]
        v[i] = s % nv
        rows[i] = offsets[i] + v[i] + 1
        s = s ÷ nv
    end
    @inbounds for i in 1:ndesc
        r = rows[i]
        @simd for j in 1:ndim
            ib[j] += cimT[j, r]
        end
    end

    @inbounds for sig in first_sig:last_sig
        # ── Fixed-point check: per-descriptor early-exit compare over ib ──
        fixed = true
        for i in 1:ndesc
            off = offsets[i]
            curm = ib[off + v[i] + 1] + margin   # unseated only by a variant beating this
            for j in 1:nvariants[i]
                if ib[off + j] > curm   # strict >: ties / within-margin keep the current variant
                    fixed = false
                    break
                end
            end
            fixed || break
        end
        fixed && push!(out, copy(v))

        # ── Odometer increment with fused row-delta ib update ──
        if sig < last_sig
            for i in 1:ndesc
                nv = nvariants[i]
                nv == 1 && continue    # radix-1: value stays 0, carry onward
                rold = rows[i]
                if v[i] + 1 < nv
                    v[i] += 1
                    rnew = rold + 1
                    rows[i] = rnew
                    @simd for j in 1:ndim
                        ib[j] += cimT[j, rnew] - cimT[j, rold]
                    end
                    break
                end
                v[i] = 0               # roll over; carry to next digit
                rnew = offsets[i] + 1
                rows[i] = rnew
                @simd for j in 1:ndim
                    ib[j] += cimT[j, rnew] - cimT[j, rold]
                end
            end
        end
    end
    return out
end

# ─── Branch-and-bound exhaustive search ─────────────────────────────────────

"""
    _bnb_bounds(cib) -> (sufmin, sufmax)

Suffix score bounds for branch-and-bound. `sufmin[k, c]` / `sufmax[k, c]`
is the minimum / maximum total contribution descriptors `k..ndesc` can make
to the impact score of flat variant column `c` over all choices of their
variants. Row `ndesc + 1` is zero, so once descriptors `1..k` are assigned,
`sufmin[k+1, c]`..`sufmax[k+1, c]` brackets what the still-free descriptors
can add to column `c`.
"""
function _bnb_bounds(cib::CIB)
    ndesc, ndim = cib.ndesc, cib.ndim
    cim_t = cib.cim_t
    sufmin = zeros(Int, ndesc + 1, ndim)
    sufmax = zeros(Int, ndesc + 1, ndim)
    @inbounds for k in ndesc:-1:1
        off = cib.desc_offsets[k]
        nv = cib.nvariants[k]
        for c in 1:ndim
            mn = typemax(Int)
            mx = typemin(Int)
            for s in 0:nv-1
                val = Int(cim_t[c, off + s + 1])   # = cim[row of variant s, c]
                mn = ifelse(val < mn, val, mn)
                mx = ifelse(val > mx, val, mx)
            end
            sufmin[k, c] = sufmin[k + 1, c] + mn
            sufmax[k, c] = sufmax[k + 1, c] + mx
        end
    end
    return sufmin, sufmax
end

# Per-task DFS state. `p` is the running prefix impact balance (exact
# contribution of the assigned descriptors to every column); `v` the current
# partial assignment. `total`/`abort`/`budget` are shared across tasks.
struct _BnBState
    cim_t::Matrix{Int}
    sufmin::Matrix{Int}
    sufmax::Matrix{Int}
    nvariants::Vector{Int}
    offsets::Vector{Int}
    ndesc::Int
    ndim::Int
    v::Vector{Int}
    p::Vector{Int}
    out::Vector{Vector{Int}}
    nodes::Base.RefValue{Int}
    total::Threads.Atomic{Int}
    abort::Threads.Atomic{Bool}
    budget::Int
    margin::Int
end

"""
    _bnb_pruned(st, k) -> Bool

With descriptors `1..k` assigned, decide whether the current subtree can be
discarded: true iff some assigned descriptor `j` has a competitor variant
whose *worst-case* completed score still exceeds the chosen variant's
*best-case* completed score by more than the rule's `margin` — then every
completion fails the fixed-point test at `j`. Ties (and gaps within the
margin) never prune, matching the favour-current convention, and at
`k == ndesc` the suffix bounds are zero, so this condition IS the exact
margin fixed-point predicate on the full scenario (`margin = 0` recovers
plain global consistency).
"""
function _bnb_pruned(st::_BnBState, k::Int)
    p = st.p
    sufmin = st.sufmin
    sufmax = st.sufmax
    kk = k + 1
    @inbounds for j in 1:k
        off = st.offsets[j]
        cs = off + st.v[j] + 1
        best_cur = p[cs] + sufmax[kk, cs] + st.margin
        for c in off+1:off+st.nvariants[j]
            if p[c] + sufmin[kk, c] > best_cur
                return true
            end
        end
    end
    return false
end

# Charge one visited node against the shared budget (flushed every 256).
# Returns false when the budget has tripped and the search must abort.
function _bnb_charge!(st::_BnBState)
    st.nodes[] += 1
    if st.nodes[] >= 256
        Threads.atomic_add!(st.total, st.nodes[])
        st.nodes[] = 0
        if st.total[] > st.budget || st.abort[]
            st.abort[] = true
            return false
        end
    end
    return true
end

# Depth-first over variants of descriptor k. Returns false on abort.
function _bnb_node!(st::_BnBState, k::Int)
    nv = st.nvariants[k]
    off = st.offsets[k]
    p = st.p
    cim_t = st.cim_t
    ndim = st.ndim
    @inbounds for s in 0:nv-1
        col = off + s + 1
        @simd for j in 1:ndim          # push chosen row onto the prefix
            p[j] += cim_t[j, col]
        end
        st.v[k] = s
        ok = _bnb_charge!(st)
        if ok && !_bnb_pruned(st, k)
            if k == st.ndesc
                push!(st.out, copy(st.v))
            else
                ok = _bnb_node!(st, k + 1)
            end
        end
        @simd for j in 1:ndim          # pop
            p[j] -= cim_t[j, col]
        end
        ok || return false
    end
    return true
end

"""
    _bnb_fixed_points(cib, sufmin, sufmax; node_budget)
        -> (Union{Nothing, Vector{Vector{Int}}}, nodes_visited)

Threaded branch-and-bound search for all fixed points. Subtrees are fanned
out as tasks over the assignments of the first `L` descriptors (chosen so
the task count comfortably exceeds the thread count). The first element is
`nothing` if the visited-node budget trips (pruning too weak to beat the
sweep), else the complete kernel sorted by ascending signature; the second
is the number of tree nodes visited (partial scenarios expanded).
"""
function _bnb_fixed_points(cib::CIB, sufmin::Matrix{Int}, sufmax::Matrix{Int};
                           node_budget::Int, margin::Int=0)
    ndesc = cib.ndesc
    nvariants = cib.nvariants
    offsets = cib.desc_offsets
    ndim = cib.ndim
    cim_t = cib.cim_t

    # Prefix depth: enough top-level assignments to keep every thread busy.
    L = 0
    nprefix = 1
    while L < ndesc && nprefix < 4 * Threads.nthreads()
        L += 1
        nprefix *= nvariants[L]
    end

    total = Threads.Atomic{Int}(0)
    abort = Threads.Atomic{Bool}(false)
    outs = [Vector{Vector{Int}}() for _ in 1:nprefix]

    @sync for pid in 0:nprefix-1
        out = outs[pid + 1]
        Threads.@spawn begin
            st = _BnBState(cim_t, sufmin, sufmax, nvariants, offsets,
                           ndesc, ndim, zeros(Int, ndesc), zeros(Int, ndim),
                           out, Ref(0), total, abort, node_budget, margin)
            # Build the prefix, checking the prune at every level so a
            # subtree dead at depth i < L is skipped without descending.
            s = pid
            alive = true
            @inbounds for i in 1:L
                st.v[i] = s % nvariants[i]
                s = s ÷ nvariants[i]
                col = offsets[i] + st.v[i] + 1
                @simd for j in 1:ndim
                    st.p[j] += cim_t[j, col]
                end
                alive = _bnb_charge!(st) && !_bnb_pruned(st, i)
                alive || break
            end
            if alive
                if L == ndesc
                    push!(out, copy(st.v))   # surviving depth-ndesc prune == fixed point
                else
                    _bnb_node!(st, L + 1)
                end
            end
            Threads.atomic_add!(total, st.nodes[])   # flush residual node count
        end
    end

    abort[] && return (nothing, total[])
    kern = Vector{Vector{Int}}()
    for out in outs
        append!(kern, out)
    end
    sort!(kern; by = u -> signature(cib, u))
    return (kern, total[])
end

"""
    find_basins(cib; rule=GlobalSuccession()) -> (fixed_points, basin_sizes, cycle_count)

Exhaustive basin-of-attraction analysis under `rule` (default
[`GlobalSuccession`](@ref)). Follows the succession chain from every scenario
in the space, counting how many starting points converge to each fixed point.

For [`GlobalSuccession`](@ref) this runs in two phases: a threaded sweep fills
a flat successor table (every scenario's succession step, computed once via an
incrementally maintained impact balance), then a sequential resolution pass
walks the table with path compression so each scenario is resolved exactly
once. Memory is ~8n bytes for the two flat tables (Int32 entries when the
space fits, Int64 otherwise), independent of the thread count. Any other rule
fills the same successor table one step at a time (still threaded, disjoint
writes) and reuses the identical resolver — so the result semantics are
identical to the fast path, only the table fill differs.

Returns:
- `fixed_points`: Vector of fixed-point scenarios (0-based variant indices),
  ordered by ascending signature
- `basin_sizes`:  corresponding basin sizes (same order as `fixed_points`)
- `cycle_count`:  number of scenarios that fall into non-fixed-point cycles
  (including scenarios whose chain merely leads into a cycle)
"""
function find_basins(cib::CIB; rule::SuccessionRule=GlobalSuccession())
    return _basins(rule, cib)
end

function _basins(::GlobalSuccession, cib::CIB)
    n = max_signature(cib) + 1
    T = _score_type(cib)
    if n <= Int(typemax(Int32)) - 1
        return _find_basins(cib, Int32, Matrix{T}(cib.cim_t))
    end
    return _find_basins(cib, Int64, Matrix{T}(cib.cim_t))
end

# Generic fallback: works for any rule. The successor of a scenario depends only
# on that scenario (the dynamics is memoryless), so the successor table can be
# filled one step at a time with disjoint, threaded writes — then handed to the
# very same path-compressed resolver the fast GlobalSuccession path uses. Only
# the table *fill* is rule-specific; resolution is shared, so semantics (into-
# cycle starts counted, not assigned to any basin) and ascending-signature
# ordering are identical to the fast path by construction. A rule may still
# override `_basins` with a faster fill of its own, but never needs to.
function _basins(rule::SuccessionRule, cib::CIB)
    n = max_signature(cib) + 1
    if n <= Int(typemax(Int32)) - 1
        return _generic_basins(rule, cib, Int32)
    end
    return _generic_basins(rule, cib, Int64)
end

function _generic_basins(rule::SuccessionRule, cib::CIB,
                         ::Type{S}) where {S<:Union{Int32,Int64}}
    n = max_signature(cib) + 1
    succ = Vector{S}(undef, n)
    nt = Threads.nthreads()
    chunk = cld(n, nt)
    @sync for t in 1:nt
        lo = (t - 1) * chunk
        hi = min(t * chunk, n) - 1
        lo > hi && continue
        Threads.@spawn for sig in lo:hi          # disjoint ranges → race-free
            u = inv_signature(cib, sig)
            v = succession_step(rule, cib, u)
            @inbounds succ[sig + 1] = S(signature(cib, v))
        end
    end
    fp_sigs, sizes, cycle_count = _resolve_and_tally(succ, n)
    fixed_points = [inv_signature(cib, s) for s in fp_sigs]
    return (fixed_points, sizes, cycle_count)
end

function _find_basins(cib::CIB, ::Type{S}, cimT::Matrix{T}) where {S<:Union{Int32,Int64}, T<:Signed}
    n = max_signature(cib) + 1
    succ = Vector{S}(undef, n)
    _successor_table!(succ, cib, cimT)
    fp_sigs, sizes, cycle_count = _resolve_and_tally(succ, n)
    kern = [inv_signature(cib, fp) for fp in fp_sigs]
    return kern, sizes, cycle_count
end

"""
    _successor_table!(succ, cib, cimT) -> succ

Fill `succ[sig + 1]` with the succession-step signature of every scenario.
Threaded over contiguous chunks (disjoint writes); each chunk advances a
mixed-radix odometer and maintains the impact balance incrementally, then
takes the per-descriptor argmax (ties favor the current variant, then the
lowest index — identical to [`succession_step`](@ref)).
"""
function _successor_table!(succ::Vector{S}, cib::CIB, cimT::Matrix{T}) where {S,T}
    n = max_signature(cib) + 1
    orders = Vector{Int}(undef, cib.ndesc)   # mixed-radix place values
    o = 1
    for i in 1:cib.ndesc
        orders[i] = o
        o *= cib.nvariants[i]
    end

    nchunks = max(1, min(n, 16 * Threads.nthreads()))
    chunk_size = cld(n, nchunks)
    nchunks = cld(n, chunk_size)
    if S === Int32 && T === Int16 && cib.ndim > 0
        # SIMD fast path: variant-major tiled argmax (identical output).
        lay = _vm_layout(cib.nvariants, cib.desc_offsets, orders, cimT)
        @sync for c in 1:nchunks
            first_sig = (c - 1) * chunk_size
            last_sig = min(c * chunk_size, n) - 1
            Threads.@spawn _successor_chunk_vm!(succ, lay, first_sig, last_sig,
                                                cib.nvariants, cib.desc_offsets,
                                                cib.ndesc, cib.ndim)
        end
        return succ
    end
    @sync for c in 1:nchunks
        first_sig = (c - 1) * chunk_size
        last_sig = min(c * chunk_size, n) - 1
        Threads.@spawn _successor_chunk!(succ, cimT, first_sig, last_sig,
                                         cib.nvariants, cib.desc_offsets,
                                         orders, cib.ndesc, cib.ndim)
    end
    return succ
end

function _successor_chunk!(succ::Vector{S}, cimT::Matrix{T},
                           first_sig::Int, last_sig::Int,
                           nvariants::Vector{Int}, offsets::Vector{Int},
                           orders::Vector{Int}, ndesc::Int, ndim::Int) where {S,T}
    v    = Vector{Int}(undef, ndesc)
    rows = Vector{Int}(undef, ndesc)
    ib   = zeros(T, ndim)

    s = first_sig
    @inbounds for i in 1:ndesc
        nv = nvariants[i]
        v[i] = s % nv
        rows[i] = offsets[i] + v[i] + 1
        s = s ÷ nv
    end
    @inbounds for i in 1:ndesc
        r = rows[i]
        @simd for j in 1:ndim
            ib[j] += cimT[j, r]
        end
    end

    @inbounds for sig in first_sig:last_sig
        # ── Successor signature: full per-descriptor argmax over ib ──
        w_sig = 0
        for i in 1:ndesc
            off = offsets[i]
            wi = v[i]
            max_val = ib[off + wi + 1]   # current variant seeds the max
            for j in 0:nvariants[i]-1
                score = ib[off + j + 1]
                better = score > max_val         # strict >: ties keep current/lower index
                max_val = ifelse(better, score, max_val)
                wi = ifelse(better, j, wi)        # branchless select (no misprediction)
            end
            w_sig += orders[i] * wi
        end
        succ[sig + 1] = S(w_sig)

        # ── Odometer increment with fused row-delta ib update ──
        if sig < last_sig
            for i in 1:ndesc
                nv = nvariants[i]
                nv == 1 && continue      # radix-1: value stays 0, carry onward
                rold = rows[i]
                if v[i] + 1 < nv
                    v[i] += 1
                    rnew = rold + 1
                    rows[i] = rnew
                    @simd for j in 1:ndim
                        ib[j] += cimT[j, rnew] - cimT[j, rold]
                    end
                    break
                end
                v[i] = 0                 # roll over; carry to next digit
                rnew = offsets[i] + 1
                rows[i] = rnew
                @simd for j in 1:ndim
                    ib[j] += cimT[j, rnew] - cimT[j, rold]
                end
            end
        end
    end
    return succ
end

# ─── SIMD successor chunk (Int16 scores, Int32 signatures) ──────────────────
#
# The scalar argmax above is the table build's hot spot (~60% of its time): a
# data-dependent scan of 3-4 variants per descriptor that cannot vectorize.
# This path restructures the impact-balance vector into *variant-major* order —
# grouping descriptors by radix so that "variant j of every descriptor in the
# group" is one contiguous plane — and then argmaxes 16 descriptors at once
# with 256-bit integer SIMD. Semantics are identical to the scalar path (the
# per-state successor table is byte-for-byte the same); only the evaluation
# order changes.

"""
    _vm_layout(nvariants, offsets, orders, cimT) -> NamedTuple

Build the variant-major layout for [`_successor_chunk_vm!`](@ref). Descriptors
are grouped by radix; within a group the balance slots are reordered so plane
`j` holds "variant `j` of each member descriptor" contiguously, and each group
is split into 16-lane tiles. Returns:

- `cimt_vm` — `cimT` with its slot axis permuted to variant-major (columns,
  indexed by the odometer's descriptor-major `rows`, are untouched);
- per-tile arrays `tile_r/tile_m/tile_base/tile_k0/tile_cur0` (radix, group
  size = plane stride, group base slot, lane offset, current-variant buffer
  offset);
- `ordbuf` — 16 `Int32` mixed-radix place values per tile, zero in pad lanes
  (pad lanes therefore contribute nothing to the signature);
- `curpos_of` — where each descriptor's current variant lives in the per-worker
  current-variant buffer (`ncur` entries).
"""
function _vm_layout(nvariants::Vector{Int}, offsets::Vector{Int},
                    orders::Vector{Int}, cimT::Matrix{Int16})
    ndesc = length(nvariants)
    ndim = size(cimT, 1)
    radixes = sort!(unique(nvariants))
    dm_of_vm = Vector{Int}(undef, ndim)      # variant-major slot -> descriptor-major slot
    curpos_of = zeros(Int, ndesc)
    tile_r = Int[]; tile_m = Int[]; tile_base = Int[]; tile_k0 = Int[]; tile_cur0 = Int[]
    ordbuf = Int32[]
    base = 0
    for r in radixes
        members = [i for i in 1:ndesc if nvariants[i] == r]
        m = length(members)
        for j in 0:r-1, (k, i) in enumerate(members)
            dm_of_vm[base + j*m + k] = offsets[i] + j + 1
        end
        for k0 in 0:16:m-1
            lanes = min(16, m - k0)
            push!(tile_r, r); push!(tile_m, m); push!(tile_base, base); push!(tile_k0, k0)
            push!(tile_cur0, length(tile_r) * 16 - 15)     # 16 cur lanes per tile
            for l in 1:16
                push!(ordbuf, l <= lanes ? Int32(orders[members[k0 + l]]) : Int32(0))
            end
            for l in 1:lanes
                curpos_of[members[k0 + l]] = (length(tile_r) - 1) * 16 + l
            end
        end
        base += r * m
    end
    return (cimt_vm = cimT[dm_of_vm, :], ntiles = length(tile_r),
            tile_r = tile_r, tile_m = tile_m, tile_base = tile_base,
            tile_k0 = tile_k0, tile_cur0 = tile_cur0, ordbuf = ordbuf,
            curpos_of = curpos_of, ncur = length(tile_r) * 16)
end

"""
    _successor_chunk_vm!(succ, lay, first_sig, last_sig, nvariants, offsets,
                         ndesc, ndim)

Variant-major successor chunk: identical output to [`_successor_chunk!`](@ref)
(same incremental odometer, same strict-`>`/current-wins tie-break), with the
per-descriptor argmax vectorized across each 16-lane tile. Per variant plane
the winner rule is a single blend:

    wins = (score > max) | (lane's current variant == j  &  score == max)

which reproduces the scalar tie-break exactly: the current variant beats an
equal incumbent (the `==` arm fires only on the lane's own current plane), all
other variants need strict `>`, and lower indices win otherwise because planes
are scanned in ascending order. Pad lanes are harmless: their loads stay within
the 16-slot padding of `ibv` and their place values in `ordbuf` are zero. The
signature is then one SIMD widening multiply + horizontal sum per tile.
"""
function _successor_chunk_vm!(succ::Vector{Int32}, lay, first_sig::Int, last_sig::Int,
                              nvariants::Vector{Int}, offsets::Vector{Int},
                              ndesc::Int, ndim::Int)
    cimt_vm = lay.cimt_vm
    v    = Vector{Int}(undef, ndesc)
    rows = Vector{Int}(undef, ndesc)
    ibv  = zeros(Int16, ndim + 16)           # +16: tile loads may overhang the last group
    cur  = zeros(Int16, lay.ncur)

    s = first_sig
    @inbounds for i in 1:ndesc
        nv = nvariants[i]
        v[i] = s % nv
        rows[i] = offsets[i] + v[i] + 1
        cur[lay.curpos_of[i]] = Int16(v[i])
        s = s ÷ nv
    end
    @inbounds for i in 1:ndesc
        r = rows[i]
        @simd for j in 1:ndim
            ibv[j] += cimt_vm[j, r]
        end
    end

    V = Vec{16,Int16}
    W = Vec{16,Int32}
    @inbounds for sig in first_sig:last_sig
        # ── Successor signature: tiled SIMD argmax over the variant planes ──
        w_sig = 0
        for t in 1:lay.ntiles
            r = lay.tile_r[t]; m = lay.tile_m[t]
            off0 = lay.tile_base[t] + lay.tile_k0[t] + 1
            cv = vload(V, cur, lay.tile_cur0[t])
            mx = vload(V, ibv, off0)                     # plane 0 seeds (lowest index)
            best = V(Int16(0))
            for j in 1:r-1
                pj = vload(V, ibv, off0 + j*m)
                jv = V(Int16(j))
                b = (pj > mx) | ((cv == jv) & (pj == mx))
                mx = vifelse(b, pj, mx)
                best = vifelse(b, jv, best)
            end
            w_sig += Int(sum(convert(W, best) * vload(W, lay.ordbuf, (t-1)*16 + 1)))
        end
        succ[sig + 1] = Int32(w_sig)

        # ── Odometer increment with fused row-delta ibv update ──
        if sig < last_sig
            for i in 1:ndesc
                nv = nvariants[i]
                nv == 1 && continue      # radix-1: value stays 0, carry onward
                rold = rows[i]
                if v[i] + 1 < nv
                    vi = v[i] + 1
                    v[i] = vi
                    rnew = rold + 1
                    rows[i] = rnew
                    cur[lay.curpos_of[i]] = Int16(vi)
                    @simd for j in 1:ndim
                        ibv[j] += cimt_vm[j, rnew] - cimt_vm[j, rold]
                    end
                    break
                end
                v[i] = 0                 # roll over; carry to next digit
                rnew = offsets[i] + 1
                rows[i] = rnew
                cur[lay.curpos_of[i]] = Int16(0)
                @simd for j in 1:ndim
                    ibv[j] += cimt_vm[j, rnew] - cimt_vm[j, rold]
                end
            end
        end
    end
    return succ
end

"""
    _fp_id!(reg_lock, fp_sig_by_id, id_by_fp_sig, sig) -> id

Locked get-or-assign of a dense 1-based id for the fixed point with signature
`sig`. Called once per fixed point discovered (≈ `nfp` times total across all
workers), so the lock is essentially uncontended. Concurrent discoverers of the
same fixed point serialize here and receive the same id.
"""
@noinline function _fp_id!(reg_lock::ReentrantLock, fp_sig_by_id::Vector{Int},
                           id_by_fp_sig::Dict{Int,Int}, sig::Int)
    lock(reg_lock)
    try
        id = get(id_by_fp_sig, sig, 0)
        if id == 0
            push!(fp_sig_by_id, sig)
            id = length(fp_sig_by_id)
            id_by_fp_sig[sig] = id
        end
        return id
    finally
        unlock(reg_lock)
    end
end

"""
    _resolve_chunk!(res, succ, lo, hi, reg_lock, fp_sig_by_id, id_by_fp_sig)

Resolve the starts in `lo:hi` into the shared `res`, following successor chains.
Cycle detection is **thread-private** (a per-worker `history` + backward scan),
so `res` only ever holds *final* labels — 0 = unvisited, -1 = cycle, k > 0 =
converges to the fixed point with dense id `k`. There is no shared in-progress
marker, so a worker that walks into another worker's not-yet-resolved chain just
re-walks it (redundant, never a false cycle) and reaches the same attractor;
every state's attractor is deterministic, so concurrent writes to the same slot
store the same value — a benign race. Fixed-point ids come from the locked
registry ([`_fp_id!`](@ref)), hit ≈ `nfp` times.
"""
function _resolve_chunk!(res::Vector{S}, succ::Vector{S}, lo::Int, hi::Int,
                         reg_lock::ReentrantLock, fp_sig_by_id::Vector{Int},
                         id_by_fp_sig::Dict{Int,Int}) where {S}
    history = Int[]
    @inbounds for start in lo:hi
        res[start + 1] != 0 && continue
        empty!(history)
        cur = start
        label = zero(S)
        while true
            r = res[cur + 1]
            if r != 0                        # already resolved: inherit (fp id or -1)
                label = r
                break
            end
            onchain = false                  # already on our own chain? -> ≥2-cycle
            for k in length(history):-1:1
                if history[k] == cur
                    onchain = true
                    break
                end
            end
            if onchain
                label = S(-1)                # whole chain (incl. pre-cycle tail) is cycle
                break
            end
            push!(history, cur)
            nxt = Int(succ[cur + 1])
            if nxt == cur                    # fixed point (counts itself: it's in history)
                label = S(_fp_id!(reg_lock, fp_sig_by_id, id_by_fp_sig, cur))
                break
            end
            cur = nxt
        end
        for h in history
            res[h + 1] = label
        end
    end
    return nothing
end

"""
    _resolve_and_tally(succ, n) -> (fp_sigs, sizes, cycle_count)

Resolve every scenario to its attractor by walking the successor table, then
tally basin sizes and the cycle count. The walk is threaded over disjoint start
ranges ([`_resolve_chunk!`](@ref)); it is race-safe because `res` holds only
final labels and every attractor is deterministic. Each fixed point gets a dense
id the first time it is reached (via a locked registry); the tally then indexes
a dense per-fixed-point counter (fixed points are few), so it costs an array
increment per scenario rather than a hash lookup. Output fixed points are sorted
by signature, so the result is identical at any thread count.
"""
function _resolve_and_tally(succ::Vector{S}, n::Int) where {S}
    res = zeros(S, n)
    reg_lock = ReentrantLock()
    fp_sig_by_id = Int[]             # dense id (1-based) -> fixed-point signature
    id_by_fp_sig = Dict{Int,Int}()   # fixed-point signature -> dense id

    nt = Threads.nthreads()
    chunk = cld(n, nt)
    @sync for t in 1:nt
        lo = (t - 1) * chunk
        hi = min(t * chunk, n) - 1
        Threads.@spawn _resolve_chunk!(res, succ, lo, hi, reg_lock,
                                       fp_sig_by_id, id_by_fp_sig)
    end

    # ── Threaded tally: res holds a dense fixed-point id (>0) or -1 (cycle) ──
    nfp = length(fp_sig_by_id)
    tally_chunk = cld(n, nt)
    local_counts = [zeros(Int, nfp) for _ in 1:nt]
    local_cyc = zeros(Int, nt)
    @sync for t in 1:nt
        lo = (t - 1) * tally_chunk
        hi = min(t * tally_chunk, n) - 1
        counts = local_counts[t]
        Threads.@spawn begin
            cyc = 0
            @inbounds for i in lo:hi
                c = Int(res[i + 1])
                if c == -1
                    cyc += 1
                else
                    counts[c] += 1              # c is a dense id in 1:nfp
                end
            end
            local_cyc[t] = cyc
        end
    end

    total = zeros(Int, nfp)
    for counts in local_counts
        @inbounds for k in 1:nfp
            total[k] += counts[k]
        end
    end
    perm = sortperm(fp_sig_by_id)
    return fp_sig_by_id[perm], total[perm], sum(local_cyc)
end

end # module

"""
    CrossImpactBalances

Cross-Impact Balance (CIB) scenario analysis. A Julia implementation of the
methodology in Weimer-Jehle (2006), ported from the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with a two-phase basin-of-attraction analysis, a threaded incremental
exhaustive sweep, and an exact branch-and-bound search that prunes
provably-inconsistent regions of the scenario space.

# Reference
Weimer-Jehle, W. (2006). Cross-impact balances: A system-theoretical approach
to cross-impact analysis. *Technological Forecasting and Social Change*,
73(4), 334–361.
"""
module CrossImpactBalances

export CIB, load_scw, load_solutions,
       impact_balance, succession_step,
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

# ─── Succession ─────────────────────────────────────────────────────────────

"""
    succession_step(cib, u) -> Vector{Int}

One step of global succession: for each descriptor, pick the variant with the
highest impact score. Ties favor the current variant, then lowest index.
"""
function succession_step(cib::CIB, u::Vector{Int})
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

# ─── Find consistent scenarios ──────────────────────────────────────────────

"""
    find_consistent(cib; algorithm=:auto) -> Vector{Vector{Int}}

Find every consistent scenario — every fixed point of the global-succession
operator — by searching the full scenario space with all available threads
(start Julia with `julia -t auto`).

`algorithm` selects the strategy:
- `:sweep` — enumerate every scenario with the incremental odometer sweep.
- `:bnb`   — branch-and-bound: assign descriptors depth-first and prune
  subtrees that provably contain no fixed point (see [`_bnb_bounds`](@ref)).
  Exact — returns the identical kernel — and typically visits a small
  fraction of the space on strongly-coupled matrices.
- `:auto`  (default) — the sweep for small spaces (< 10^5 scenarios),
  otherwise branch-and-bound with a node budget of `n ÷ 16`; if pruning is
  too weak to pay off, the budget trips and the sweep runs instead.

The returned kernel is ordered by ascending signature for every algorithm.
Only true fixed points are returned; scenarios that fall into longer cycles
never appear.
"""
function find_consistent(cib::CIB; algorithm::Symbol=:auto,
                         bnb_node_budget::Union{Nothing,Int}=nothing)
    algorithm in (:auto, :bnb, :sweep) ||
        throw(ArgumentError("algorithm must be :auto, :bnb or :sweep, got $(repr(algorithm))"))
    n = max_signature(cib) + 1
    if algorithm == :sweep || (algorithm == :auto && n < 100_000)
        return _find_consistent_exhaustive(cib)
    end
    sufmin, sufmax = _bnb_bounds(cib)
    budget = something(bnb_node_budget,
                       algorithm == :bnb ? typemax(Int) : n ÷ 16)
    kern, _ = _bnb_fixed_points(cib, sufmin, sufmax; node_budget=budget)
    !isnothing(kern) && return kern
    # Budget tripped: pruning too weak on this matrix — fall back.
    return _find_consistent_exhaustive(cib)
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
function _find_consistent_exhaustive(cib::CIB)
    T = _score_type(cib)
    return _sweep_fixed_points(cib, Matrix{T}(cib.cim_t))
end

function _sweep_fixed_points(cib::CIB, cimT::Matrix{T}) where {T<:Signed}
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
                                     cib.ndesc, cib.ndim)
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
                       ndesc::Int, ndim::Int) where {T<:Signed}
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
            cur = ib[off + v[i] + 1]
            for j in 1:nvariants[i]
                if ib[off + j] > cur   # strict >: ties keep the current variant
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
end

"""
    _bnb_pruned(st, k) -> Bool

With descriptors `1..k` assigned, decide whether the current subtree can be
discarded: true iff some assigned descriptor `j` has a competitor variant
whose *worst-case* completed score still strictly exceeds the chosen
variant's *best-case* completed score — then every completion fails the
fixed-point test at `j`. Ties never prune (matching the favour-current
convention), and at `k == ndesc` the suffix bounds are zero, so this
condition IS the exact fixed-point predicate on the full scenario.
"""
function _bnb_pruned(st::_BnBState, k::Int)
    p = st.p
    sufmin = st.sufmin
    sufmax = st.sufmax
    kk = k + 1
    @inbounds for j in 1:k
        off = st.offsets[j]
        cs = off + st.v[j] + 1
        best_cur = p[cs] + sufmax[kk, cs]
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
                           node_budget::Int)
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
                           out, Ref(0), total, abort, node_budget)
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
    find_basins(cib) -> (fixed_points, basin_sizes, cycle_count)

Exhaustive basin-of-attraction analysis. Follows the succession chain from
every scenario in the space, counting how many starting points converge to
each fixed point.

Runs in two phases: a threaded sweep fills a flat successor table (every
scenario's succession step, computed once via an incrementally maintained
impact balance), then a sequential resolution pass walks the table with
path compression so each scenario is resolved exactly once. Memory is
~8n bytes for the two flat tables (Int32 entries when the space fits,
Int64 otherwise), independent of the thread count.

Returns:
- `fixed_points`: Vector of fixed-point scenarios (0-based variant indices),
  ordered by ascending signature
- `basin_sizes`:  corresponding basin sizes (same order as `fixed_points`)
- `cycle_count`:  number of scenarios that fall into non-fixed-point cycles
  (including scenarios whose chain merely leads into a cycle)
"""
function find_basins(cib::CIB)
    n = max_signature(cib) + 1
    T = _score_type(cib)
    if n <= Int(typemax(Int32)) - 1
        return _find_basins(cib, Int32, Matrix{T}(cib.cim_t))
    end
    return _find_basins(cib, Int64, Matrix{T}(cib.cim_t))
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
                if score > max_val       # strict >: ties keep current/lower index
                    max_val = score
                    wi = j
                end
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

"""
    _resolve_and_tally(succ, n) -> (fp_sigs, sizes, cycle_count)

Resolve every scenario to its attractor by walking the successor table.
`res` values: 0 = unvisited, -2 = on the chain currently being walked,
-1 = falls into (or leads into) a non-fixed-point cycle, k > 0 = converges
to the fixed point with signature k - 1. The in-progress marker makes cycle
detection O(1) per step; on any resolution the entire walked chain is
back-filled, so total work is O(n). The walk is sequential (the table is
shared mutable state); the final tally pass is threaded.
"""
function _resolve_and_tally(succ::Vector{S}, n::Int) where {S}
    res = zeros(S, n)
    history = Int[]
    marker = S(-2)

    @inbounds for start in 0:n-1
        res[start + 1] != 0 && continue
        empty!(history)
        cur = start
        res[cur + 1] = marker
        push!(history, cur)

        while true
            nxt = Int(succ[cur + 1])
            r = res[nxt + 1]
            if r == marker
                # Hit our own chain: either the last element is a fixed
                # point (succ maps it to itself) or we closed a true cycle.
                # Chains *leading into* a cycle count as cycle scenarios too.
                val = nxt == cur ? S(cur + 1) : S(-1)
                for h in history
                    res[h + 1] = val
                end
                break
            elseif r != 0
                for h in history        # already resolved: backfill chain
                    res[h + 1] = r
                end
                break
            end
            res[nxt + 1] = marker
            push!(history, nxt)
            cur = nxt
        end
    end

    # ── Threaded tally (read-only over res) ──
    nt = Threads.nthreads()
    tally_chunk = cld(n, nt)
    local_counts = [Dict{Int,Int}() for _ in 1:nt]
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
                    counts[c - 1] = get(counts, c - 1, 0) + 1
                end
            end
            local_cyc[t] = cyc
        end
    end

    merged = Dict{Int,Int}()
    for counts in local_counts
        for (fp, cnt) in counts
            merged[fp] = get(merged, fp, 0) + cnt
        end
    end
    fp_sigs = sort!(collect(keys(merged)))
    return fp_sigs, [merged[fp] for fp in fp_sigs], sum(local_cyc)
end

end # module

"""
    CrossImpactBalances

Cross-Impact Balance (CIB) scenario analysis. A Julia implementation of the
methodology in Weimer-Jehle (2006), ported from the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with a memoized basin-of-attraction analysis and a threaded exhaustive search.

# Reference
Weimer-Jehle, W. (2006). Cross-impact balances: A system-theoretical approach
to cross-impact analysis. *Technological Forecasting and Social Change*,
73(4), 334–361.
"""
module CrossImpactBalances

using Random
using SparseArrays

export CIB, load_scw, load_solutions,
       impact_balance, own_impact_balance, cross_impact_balance, inner_product,
       succession_step, succession, find_consistent, find_basins,
       signature, inv_signature, max_signature,
       sim_anneal, inner_product_matrix, build_graph, merge_scenarios,
       set_thresholds!, rand_scenario

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
- `thresholds`: per-descriptor thresholds for simulated annealing
- `mc_threshold`: cutoff for switching to Monte Carlo sampling
- `desc_offsets`: cumulative 0-based variant offsets per descriptor

`CIB` is immutable in its field bindings. `kernel` and `thresholds` are
`Vector`s whose contents may be mutated in place (see `set_thresholds!`).
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
    thresholds::Vector{Int}
    mc_threshold::Int
    desc_offsets::Vector{Int}
end

# ─── .scw file parser ───────────────────────────────────────────────────────

"""
    load_scw(scw_file; sl_file=nothing, kernel=nothing,
             mc_threshold=10000, exhaustive=false,
             rng=Random.default_rng()) -> CIB

Parse a ScenarioWizard .scw file and optionally a .sl solutions file.
Returns a fully populated CIB object.

When `exhaustive=true`, every scenario in the space is checked (using threads
when Julia is started with multiple threads, e.g. `julia -t8`).

When the scenario space exceeds `mc_threshold` and `exhaustive=false`, a
Monte-Carlo sample of that size is drawn without replacement using `rng`.
"""
function load_scw(scw_file::String; sl_file::Union{String,Nothing}=nothing,
                  kernel::Union{Vector{Vector{Int}},Nothing}=nothing,
                  mc_threshold::Int=10000,
                  exhaustive::Bool=false,
                  rng::AbstractRNG=Random.default_rng())
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
    thresholds = zeros(Int, ndesc)

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
              Vector{Vector{Int}}(), thresholds, mc_threshold, desc_offsets)

    if !isnothing(kernel)
        append!(cib.kernel, kernel)
    elseif !isnothing(sl_file)
        append!(cib.kernel, load_solutions(cib, sl_file))
    else
        append!(cib.kernel, find_consistent(cib; exhaustive=exhaustive, rng=rng))
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

# ─── Index conversion helpers ───────────────────────────────────────────────

"""
    varndx_to_tablendx(cib, u) -> Vector{Int}

Convert a scenario `u` (length-`ndesc` vector of 0-based variant indices) to
the corresponding 1-based row/column indices into `cib.cim`, i.e. the indices
of `u`'s selected variants in the flat variant space.

This is an internal helper used by [`impact_balance`](@ref) and friends. Hot
loops in `find_basins` and `_find_consistent_exhaustive` inline the equivalent
arithmetic to avoid the per-call allocation here.
"""
function varndx_to_tablendx(cib::CIB, u::Vector{Int})
    tablendx = similar(u)
    offset = 0
    for i in eachindex(u)
        tablendx[i] = offset + u[i] + 1  # +1 for Julia 1-based indexing
        offset += cib.nvariants[i]
    end
    return tablendx
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

# ─── Scenario sampling ──────────────────────────────────────────────────────

"""
    get_scenario_signatures(cib; max=cib.mc_threshold, allow_dups=false, rng=...)

Return scenario signatures: every signature in `0:n-1` if the space size `n`
does not exceed `max`, otherwise a random sample of `max` signatures.

By default the sample is drawn without replacement, matching Python CIBSA's
`np.random.choice(..., replace=False)` behavior. Pass `allow_dups=true` to
sample with replacement (cheaper but biased toward duplicates).

`rng` controls reproducibility; pass a seeded `AbstractRNG` for deterministic
sampling.
"""
function get_scenario_signatures(cib::CIB;
                                 max::Int=cib.mc_threshold,
                                 allow_dups::Bool=false,
                                 rng::AbstractRNG=Random.default_rng())
    n = max_signature(cib) + 1
    if n <= max
        return 0:n-1
    end
    if allow_dups
        return rand(rng, 0:n-1, max)
    end
    return _sample_without_replacement(rng, n, max)
end

# Reservoir-style partial Fisher-Yates: O(max) extra memory + O(max) swaps.
# Returns `max` distinct integers from 0:n-1 in random order.
function _sample_without_replacement(rng::AbstractRNG, n::Int, max::Int)
    @assert max <= n
    # When `max` is close to `n` a dense permutation is cheaper; otherwise a
    # hash-table-based partial shuffle keeps memory at O(max).
    if max * 4 >= n
        perm = randperm(rng, n) .- 1  # 0-based signatures
        return perm[1:max]
    end
    seen = Dict{Int,Int}()  # virtual swap map: i -> what value lives at position i
    result = Vector{Int}(undef, max)
    @inbounds for i in 1:max
        # Pick j uniformly from {i, i+1, ..., n} (1-based positions in conceptual array)
        j = rand(rng, i:n)
        # Read the value at position j (default = j-1, the original integer)
        vj = get(seen, j, j - 1)
        # Read the value at position i (default = i-1)
        vi = get(seen, i, i - 1)
        # Swap positions i and j conceptually
        seen[j] = vi
        result[i] = vj
    end
    return result
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

"""
    own_impact_balance(cib, u) -> Vector{Int}

The impact balance entries corresponding to the scenario's own selected variants.
"""
function own_impact_balance(cib::CIB, u::Vector{Int})
    ib = impact_balance(cib, u)
    return ib[varndx_to_tablendx(cib, u)]
end

"""
    cross_impact_balance(cib, u, v) -> Vector{Int}

The impact balance of scenario `u` evaluated at the variants of scenario `v`.
"""
function cross_impact_balance(cib::CIB, u::Vector{Int}, v::Vector{Int})
    ib = impact_balance(cib, u)
    return ib[varndx_to_tablendx(cib, v)]
end

"""
    inner_product(cib, u, v) -> Int

Inner product of impact balance of `u` evaluated at `v`.
"""
function inner_product(cib::CIB, u::Vector{Int}, v::Vector{Int})
    ib = impact_balance(cib, u)
    return sum(ib[varndx_to_tablendx(cib, v)])
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

"""
    succession(cib, u) -> (cycle_length, attractor)

Follow global succession from scenario `u` until convergence to a fixed point
or detection of a cycle. Returns (cycle_length, final_scenario).
cycle_length=1 means a consistent scenario (fixed point).

Cycle detection uses an O(1)-per-step hashtable, so total work is linear in
the trajectory length.
"""
function succession(cib::CIB, u::Vector{Int})
    start_sig = signature(cib, u)
    history_sig = Int[start_sig]
    seen = Dict{Int,Int}()  # signature -> 1-based position in history_sig
    seen[start_sig] = 1
    v = copy(u)
    while true
        v = succession_step(cib, v)
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
    find_consistent(cib; ignore_cycles=true, exhaustive=false,
                    rng=Random.default_rng()) -> Vector{Vector{Int}}

Find all consistent scenarios by running succession from every (or sampled) starting point.

When `exhaustive=true`, every scenario in the space is enumerated, using all
available threads (start Julia with `julia -t auto` or `julia -t8`).
This guarantees finding all fixed points but requires iterating the full space.

When `exhaustive=false` and the space exceeds `cib.mc_threshold`, scenarios
are sampled without replacement using `rng`.
"""
function find_consistent(cib::CIB; ignore_cycles::Bool=true, exhaustive::Bool=false,
                         rng::AbstractRNG=Random.default_rng())
    if exhaustive
        return _find_consistent_exhaustive(cib; ignore_cycles=ignore_cycles)
    end
    kern = Vector{Vector{Int}}()
    seen = Set{Int}()
    for v_sig in get_scenario_signatures(cib; rng=rng)
        v = inv_signature(cib, v_sig)
        nper, veqm = succession(cib, v)
        if ignore_cycles && nper > 1
            continue
        end
        veqm_sig = signature(cib, veqm)
        if !(veqm_sig in seen)
            push!(seen, veqm_sig)
            push!(kern, veqm)
        end
    end
    return kern
end

"""
    _find_consistent_exhaustive(cib; ignore_cycles=true) -> Vector{Vector{Int}}

Internal threaded exhaustive search. Partitions the signature space across
`Threads.nthreads()` chunks; within each chunk uses a mixed-radix counter to
avoid a `divmod` per signature, and a per-descriptor early-exit fixed-point
check that bails as soon as any competitor variant beats the current one.

`ignore_cycles` is accepted for API symmetry with [`find_consistent`](@ref)
but is unused here: this routine only emits true fixed points (cycle length
1), so non-fixed-point cycles never appear.
"""
function _find_consistent_exhaustive(cib::CIB; ignore_cycles::Bool=true)
    n = max_signature(cib) + 1
    nt = Threads.nthreads()
    chunk_size = cld(n, nt)
    ndesc = cib.ndesc
    cim = cib.cim
    nvariants = cib.nvariants
    offsets = cib.desc_offsets

    local_kerns = [Vector{Vector{Int}}() for _ in 1:nt]

    Threads.@threads for chunk in 1:nt
        v    = Vector{Int}(undef, ndesc)
        tndx = Vector{Int}(undef, ndesc)

        first_sig = (chunk - 1) * chunk_size
        last_sig  = min(chunk * chunk_size - 1, n - 1)

        # Decode first_sig once (the only divmod in this chunk)
        s = first_sig
        @inbounds for i in 1:ndesc
            nv = nvariants[i]
            v[i] = s % nv
            tndx[i] = offsets[i] + v[i] + 1
            s = s ÷ nv
        end

        for sig in first_sig:last_sig
            # ── Per-descriptor fixed-point check with early exit ──
            # Compute impact score only for columns we need. If any
            # competitor variant beats the current one, bail immediately.
            fixed = true
            @inbounds for i in 1:ndesc
                nv = nvariants[i]
                off = offsets[i]
                cur_col = off + v[i] + 1
                # Score of current variant
                cur_score = 0
                for k in 1:ndesc
                    cur_score += cim[tndx[k], cur_col]
                end
                # Does any other variant beat it?
                for j in 0:nv-1
                    j == v[i] && continue
                    col = off + j + 1
                    score = 0
                    for k in 1:ndesc
                        score += cim[tndx[k], col]
                    end
                    if score > cur_score
                        fixed = false
                        break
                    end
                end
                fixed || break
            end

            if fixed
                push!(local_kerns[chunk], copy(v))
            end

            # ── Mixed-radix increment (replaces inv_signature divmod) ──
            if sig < last_sig
                @inbounds for i in 1:ndesc
                    v[i] += 1
                    tndx[i] += 1
                    if v[i] < nvariants[i]
                        break
                    end
                    v[i] = 0
                    tndx[i] = offsets[i] + 1
                end
            end
        end
    end

    # Merge chunk-local results (no dedup needed — each point is unique)
    kern = Vector{Vector{Int}}()
    for chunk in 1:nt
        append!(kern, local_kerns[chunk])
    end
    return kern
end

"""
    find_basins(cib) -> (fixed_points, basin_sizes, cycle_count)

Exhaustive basin-of-attraction analysis. Follows the succession chain from
every scenario in the space, counting how many starting points converge to
each fixed point.

Uses a flat cache so each scenario is resolved exactly once (O(n) total
succession steps). Returns:
- `fixed_points`: Vector of fixed-point scenarios (0-based variant indices)
- `basin_sizes`:  corresponding basin sizes (same order as `fixed_points`)
- `cycle_count`:  number of scenarios that fall into non-fixed-point cycles
"""
function find_basins(cib::CIB)
    n = max_signature(cib) + 1
    nt = Threads.nthreads()
    ndesc = cib.ndesc
    ndim = cib.ndim
    cim_t = cib.cim_t   # row-vectors of cim live as columns of cim_t (SIMD-friendly)
    nvariants = cib.nvariants
    offsets = cib.desc_offsets

    chunk_size = cld(n, nt)

    # Per-thread caches (memory cost: nt * n * 8 bytes — matches the single-thread
    # version's footprint times the thread count). Each thread also produces a
    # local basin tally; we merge at the end.
    # cache values: 0 = unvisited, -1 = cycle, k>0 = converges to fp with sig (k-1).
    thread_basins = [Dict{Int,Int}() for _ in 1:nt]
    thread_cycle_counts = zeros(Int, nt)

    Threads.@threads for tid in 1:nt
        first_sig = (tid - 1) * chunk_size
        last_sig  = min(tid * chunk_size, n) - 1
        first_sig > last_sig && continue

        cache    = zeros(Int, n)
        w        = Vector{Int}(undef, ndesc)
        ib       = Vector{Int}(undef, ndim)
        tndx     = Vector{Int}(undef, ndesc)
        history  = Int[]

        for start_sig in first_sig:last_sig
            @inbounds cache[start_sig + 1] != 0 && continue

            # Decode start_sig into w
            s = start_sig
            @inbounds for i in 1:ndesc
                w[i] = s % nvariants[i]
                s = s ÷ nvariants[i]
            end

            empty!(history)
            push!(history, start_sig)

            while true
                # ── inline succession_step: row-at-a-time with @simd ──
                @inbounds for i in 1:ndesc
                    tndx[i] = offsets[i] + w[i] + 1
                end
                r1 = tndx[1]
                @inbounds @simd for j in 1:ndim
                    ib[j] = cim_t[j, r1]
                end
                for ki in 2:ndesc
                    r = @inbounds tndx[ki]
                    @inbounds @simd for j in 1:ndim
                        ib[j] += cim_t[j, r]
                    end
                end
                # Pick best variant per descriptor
                @inbounds for i in 1:ndesc
                    nv = nvariants[i]
                    off = offsets[i]
                    max_val = ib[off + w[i] + 1]
                    for j in 0:nv-1
                        score = ib[off + j + 1]
                        if score > max_val
                            max_val = score
                            w[i] = j
                        end
                    end
                end

                # Signature of w (the successor)
                w_sig = 0
                order = 1
                @inbounds for i in 1:ndesc
                    w_sig += order * w[i]
                    order *= nvariants[i]
                end

                # ── Already resolved? Backfill entire history ──
                @inbounds cached = cache[w_sig + 1]
                if cached != 0
                    @inbounds for h in history
                        cache[h + 1] = cached
                    end
                    break
                end

                # ── Cycle detection: is w_sig already in our chain? ──
                cycle_start = 0
                for k in length(history):-1:1
                    if @inbounds history[k] == w_sig
                        cycle_start = k
                        break
                    end
                end

                if cycle_start > 0
                    if cycle_start == length(history)
                        val = w_sig + 1
                        @inbounds for h in history
                            cache[h + 1] = val
                        end
                    else
                        @inbounds for h in history
                            cache[h + 1] = -1
                        end
                    end
                    break
                end

                push!(history, w_sig)
            end
        end

        # ── Tally THIS thread's chunk through its own cache ──
        local_basins = thread_basins[tid]
        local_cycle = 0
        @inbounds for sig in first_sig:last_sig
            c = cache[sig + 1]
            if c == -1
                local_cycle += 1
            elseif c > 0
                fp_sig = c - 1
                local_basins[fp_sig] = get(local_basins, fp_sig, 0) + 1
            end
        end
        thread_cycle_counts[tid] = local_cycle
    end

    # ── Merge per-thread basin tallies ──
    cycle_count = sum(thread_cycle_counts)
    kern = Vector{Vector{Int}}()
    basins = Vector{Int}()
    fp_order = Dict{Int, Int}()
    for tid in 1:nt
        for (fp_sig, count) in thread_basins[tid]
            idx = get(fp_order, fp_sig, 0)
            if idx == 0
                push!(kern, inv_signature(cib, fp_sig))
                push!(basins, count)
                fp_order[fp_sig] = length(kern)
            else
                @inbounds basins[idx] += count
            end
        end
    end

    return kern, basins, cycle_count
end

# ─── Simulated annealing ────────────────────────────────────────────────────

"""
    sim_anneal(cib, u; ignore_cycles=true, return_weights=false,
               rng=Random.default_rng())

Threshold-gated accessibility analysis (a misnomer carried over from the
Python CIBSA reference — it is *not* a classical simulated-annealing
optimizer). For each candidate scenario `v` drawn from the scenario space (or
a Monte-Carlo sample when it exceeds `cib.mc_threshold`), accept `v` only if
its cross-impact at `v` exceeds `u`'s own impact by at least
`cib.thresholds[i]` in every descriptor `i`. Each accepted `v` is then
deterministically run to a fixed point via [`succession`](@ref).

When `return_weights=true`, returns a `Dict` mapping the signature of each
reachable fixed point to the number of accepted candidates that converge to
it (with a special `"reject"` key for candidates that fail the threshold
gate). Otherwise returns the list of distinct reachable fixed points.

`rng` controls the Monte-Carlo sample when one is drawn.
"""
function sim_anneal(cib::CIB, u::Vector{Int};
                    ignore_cycles::Bool=true, return_weights::Bool=false,
                    rng::AbstractRNG=Random.default_rng())
    accessible = Vector{Vector{Int}}()
    weights = Dict{Any,Int}()
    weights["reject"] = 0
    uib = own_impact_balance(cib, u)
    sig_set = get_scenario_signatures(cib; rng=rng)

    for v_sig in sig_set
        v = inv_signature(cib, v_sig)
        xib = cross_impact_balance(cib, u, v)
        valid = true
        for (ui, xi, thr) in zip(uib, xib, cib.thresholds)
            if xi + thr <= ui
                valid = false
                break
            end
        end
        if valid
            nper, veqm = succession(cib, v)
            if ignore_cycles && nper > 1
                continue
            end
            veqm_sig = signature(cib, veqm)
            if haskey(weights, veqm_sig)
                weights[veqm_sig] += 1
            else
                weights[veqm_sig] = 1
                push!(accessible, veqm)
            end
        else
            weights["reject"] += 1
        end
    end

    return return_weights ? weights : accessible
end

# ─── Inner product matrix ───────────────────────────────────────────────────

"""
    inner_product_matrix(cib) -> Matrix{Int}

Compute the matrix of inner products between all kernel scenarios.

Threaded over rows when Julia is started with multiple threads. Each row
hoists its single [`impact_balance`](@ref) call out of the inner loop, so
total work is `O(k × ndim)` impact-balance accumulations rather than `O(k²)`.
"""
function inner_product_matrix(cib::CIB)
    k = length(cib.kernel)
    M = zeros(Int, k, k)
    # Precompute table indices for every kernel element once
    tndxs = [varndx_to_tablendx(cib, v) for v in cib.kernel]
    Threads.@threads for i in 1:k
        ib = impact_balance(cib, cib.kernel[i])
        @inbounds for j in 1:k
            s = 0
            for r in tndxs[j]
                s += ib[r]
            end
            M[i, j] = s
        end
    end
    return M
end

# ─── Graph and merge ────────────────────────────────────────────────────────

"""
    build_graph(cib) -> SparseMatrixCSC

Build an adjacency matrix of kernel scenarios connected through simulated annealing.
"""
function build_graph(cib::CIB)
    k = length(cib.kernel)
    kernel_sigs = [signature(cib, u) for u in cib.kernel]
    adj = zeros(Int, k, k)
    for (r, u) in enumerate(cib.kernel)
        accessible = sim_anneal(cib, u)
        for w in accessible
            w_sig = signature(cib, w)
            idx = findfirst(==(w_sig), kernel_sigs)
            if !isnothing(idx)
                adj[r, idx] = 1
            end
        end
    end
    return sparse(adj)
end

"""
    merge_scenarios(cib) -> Vector{Vector{Int}}

Find connected components of the fluctuation-connected graph.
Returns groups of scenario signatures that merge under simulated annealing.
"""
function merge_scenarios(cib::CIB)
    adj = build_graph(cib)
    k = size(adj, 1)

    # Simple BFS connected components (avoids needing Graphs.jl)
    visited = falses(k)
    components = Vector{Vector{Int}}()

    for start in 1:k
        visited[start] && continue
        component = Int[]
        queue = [start]
        visited[start] = true
        while !isempty(queue)
            node = popfirst!(queue)
            push!(component, signature(cib, cib.kernel[node]))
            for neighbor in 1:k
                if !visited[neighbor] && (adj[node, neighbor] != 0 || adj[neighbor, node] != 0)
                    visited[neighbor] = true
                    push!(queue, neighbor)
                end
            end
        end
        push!(components, component)
    end

    return components
end

# ─── Mutators and helpers ────────────────────────────────────────────────────

"""
    set_thresholds!(cib::CIB, thr::AbstractVector{<:Integer}) -> CIB

Replace `cib.thresholds` element-wise with `thr`. The length of `thr` must
equal `cib.ndesc`. Returns `cib` so calls can be chained.

The threshold vector is consumed by [`sim_anneal`](@ref), [`build_graph`](@ref),
and [`merge_scenarios`](@ref); it has no effect on `find_consistent` or
`find_basins`.
"""
function set_thresholds!(cib::CIB, thr::AbstractVector{<:Integer})
    length(thr) == cib.ndesc ||
        throw(DimensionMismatch("set_thresholds!: expected length $(cib.ndesc), got $(length(thr))"))
    copyto!(cib.thresholds, thr)
    return cib
end

"""
    rand_scenario(cib::CIB; rng=Random.default_rng()) -> Vector{Int}

Return a uniformly random scenario, i.e. a length-`ndesc` vector of 0-based
variant indices with `0 ≤ u[i] < cib.nvariants[i]`.
"""
function rand_scenario(cib::CIB; rng::AbstractRNG=Random.default_rng())
    u = Vector{Int}(undef, cib.ndesc)
    @inbounds for i in 1:cib.ndesc
        u[i] = rand(rng, 0:cib.nvariants[i] - 1)
    end
    return u
end

end # module

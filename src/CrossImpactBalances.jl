module CrossImpactBalances

using SparseArrays

export CIB, load_scw, load_solutions,
       impact_balance, own_impact_balance, cross_impact_balance, inner_product,
       succession_step, succession, find_consistent,
       signature, inv_signature, max_signature,
       sim_anneal, inner_product_matrix, build_graph, merge_scenarios

"""
    CIB

Cross-impact balance analysis object.

Fields:
- `descriptors`: ordered list of descriptor names
- `variants`: dict mapping descriptor name to list of variant names
- `nvariants`: number of variants per descriptor
- `cim`: the cross-impact matrix (n × n), n = sum of all variants
- `ndim`: total number of variants (size of CIM)
- `ndesc`: number of descriptors
- `kernel`: list of consistent scenarios (each a Vector{Int}, 0-based variant indices)
- `thresholds`: per-descriptor thresholds for simulated annealing
- `mc_threshold`: cutoff for switching to Monte Carlo sampling
"""
mutable struct CIB
    descriptors::Vector{String}
    variants::Dict{String, Vector{String}}
    nvariants::Vector{Int}
    cim::Matrix{Int}
    ndim::Int
    ndesc::Int
    kernel::Vector{Vector{Int}}
    thresholds::Vector{Int}
    mc_threshold::Int
end

# ─── .scw file parser ───────────────────────────────────────────────────────

"""
    load_scw(scw_file; sl_file=nothing, kernel=nothing, mc_threshold=10000, exhaustive=false) -> CIB

Parse a ScenarioWizard .scw file and optionally a .sl solutions file.
Returns a fully populated CIB object.

When `exhaustive=true`, every scenario in the space is checked (using threads
when Julia is started with multiple threads, e.g. `julia -t8`).
"""
function load_scw(scw_file::String; sl_file::Union{String,Nothing}=nothing,
                  kernel::Union{Vector{Vector{Int}},Nothing}=nothing,
                  mc_threshold::Int=10000,
                  exhaustive::Bool=false)
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

    # Build CIM matrix
    cim = zeros(Int, n, n)
    for (i, row) in enumerate(cim_rows)
        for (j, val) in enumerate(row)
            cim[i, j] = val
        end
    end

    ndesc = d + 1
    thresholds = zeros(Int, ndesc)

    cib = CIB(descriptors, variants, nvars, cim, n, ndesc,
              Vector{Vector{Int}}(), thresholds, mc_threshold)

    # Load kernel
    if !isnothing(kernel)
        cib.kernel = kernel
    elseif !isnothing(sl_file)
        cib.kernel = load_solutions(cib, sl_file)
    else
        cib.kernel = find_consistent(cib; exhaustive=exhaustive)
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
Convert 0-based variant indices (length ndesc) to 1-based CIM row/column indices.
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
Return scenario signatures: all of them, or a random sample if the space exceeds mc_threshold.
"""
function get_scenario_signatures(cib::CIB)
    n = max_signature(cib) + 1
    if n > cib.mc_threshold
        return rand(0:n-1, cib.mc_threshold)
    else
        return 0:n-1
    end
end

# ─── Impact balance ─────────────────────────────────────────────────────────

"""
    impact_balance(cib, u) -> Vector{Int}

Compute the impact balance vector for scenario `u`.
Returns a vector of length ndim (one score per variant across all descriptors).
"""
function impact_balance(cib::CIB, u::Vector{Int})
    rows = varndx_to_tablendx(cib, u)
    # Sum the selected rows of the CIM
    ib = zeros(Int, cib.ndim)
    for r in rows
        for j in 1:cib.ndim
            ib[j] += cib.cim[r, j]
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
"""
function succession(cib::CIB, u::Vector{Int})
    iterations_sig = Int[signature(cib, u)]
    v = copy(u)
    while true
        v = succession_step(cib, v)
        v_sig = signature(cib, v)
        # Check history in reverse for cycle detection
        n = 1
        foundit = false
        for k in length(iterations_sig):-1:1
            if iterations_sig[k] == v_sig
                foundit = true
                break
            end
            n += 1
        end
        if foundit
            return (n, v)
        end
        push!(iterations_sig, v_sig)
    end
end

# ─── Find consistent scenarios ──────────────────────────────────────────────

"""
    find_consistent(cib; ignore_cycles=true, exhaustive=false) -> Vector{Vector{Int}}

Find all consistent scenarios by running succession from every (or sampled) starting point.

When `exhaustive=true`, every scenario in the space is enumerated, using all
available threads (start Julia with `julia -t auto` or `julia -t8`).
This guarantees finding all fixed points but requires iterating the full space.
"""
function find_consistent(cib::CIB; ignore_cycles::Bool=true, exhaustive::Bool=false)
    if exhaustive
        return _find_consistent_exhaustive(cib; ignore_cycles=ignore_cycles)
    end
    kern = Vector{Vector{Int}}()
    seen = Set{Int}()
    for v_sig in get_scenario_signatures(cib)
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

function _find_consistent_exhaustive(cib::CIB; ignore_cycles::Bool=true)
    n = max_signature(cib) + 1
    nt = Threads.nthreads()
    chunk_size = cld(n, nt)
    ndesc = cib.ndesc
    ndim = cib.ndim

    local_kerns = [Vector{Vector{Int}}() for _ in 1:nt]

    Threads.@threads for chunk in 1:nt
        # Pre-allocated working buffers — zero allocations in the hot loop
        v  = Vector{Int}(undef, ndesc)
        ib = Vector{Int}(undef, ndim)
        tndx = Vector{Int}(undef, ndesc)

        first_sig = (chunk - 1) * chunk_size
        last_sig  = min(chunk * chunk_size - 1, n - 1)

        for sig in first_sig:last_sig
            # inv_signature in-place
            s = sig
            @inbounds for i in 1:ndesc
                nv = cib.nvariants[i]
                v[i] = s % nv
                s = s ÷ nv
            end

            # One succession_step: compute impact balance, pick best variant
            offset = 0
            @inbounds for i in 1:ndesc
                tndx[i] = offset + v[i] + 1
                offset += cib.nvariants[i]
            end
            fill!(ib, 0)
            @inbounds for r in tndx
                for j in 1:ndim
                    ib[j] += cib.cim[r, j]
                end
            end

            # Check if v is a fixed point: does each descriptor already hold
            # the variant with the highest impact score?
            fixed = true
            start = 1
            @inbounds for i in 1:ndesc
                nv = cib.nvariants[i]
                max_val = ib[start + v[i]]
                for j in 0:nv-1
                    if ib[start + j] > max_val
                        fixed = false
                        break
                    end
                end
                fixed || break
                start += nv
            end

            if fixed
                push!(local_kerns[chunk], copy(v))
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

# ─── Simulated annealing ────────────────────────────────────────────────────

"""
    sim_anneal(cib, u; ignore_cycles=true, return_weights=false)

Simulated annealing: find scenarios accessible from `u` under threshold perturbation.
If return_weights=true, returns a Dict mapping signature => count.
Otherwise returns a list of accessible attractor scenarios.
"""
function sim_anneal(cib::CIB, u::Vector{Int};
                    ignore_cycles::Bool=true, return_weights::Bool=false)
    accessible = Vector{Vector{Int}}()
    weights = Dict{Any,Int}()
    weights["reject"] = 0
    uib = own_impact_balance(cib, u)
    sig_set = get_scenario_signatures(cib)

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
"""
function inner_product_matrix(cib::CIB)
    k = length(cib.kernel)
    M = zeros(Int, k, k)
    for (i, u) in enumerate(cib.kernel)
        for (j, v) in enumerate(cib.kernel)
            M[i, j] = inner_product(cib, u, v)
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

end # module

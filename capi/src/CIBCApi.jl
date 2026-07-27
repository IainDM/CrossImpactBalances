"""
    CIBCApi

A thin C-callable surface over `CrossImpactBalances`, intended to be compiled
into a standalone shared library with PackageCompiler (`create_library`). It
lets non-Julia callers (e.g. a Python `ctypes` wrapper) drive the engine
without shipping any Julia source.

Design
------
* Every exported `cib_*` function takes and returns **JSON** encoded as C
  strings (`Cstring`). Returned strings are heap-allocated with `Libc.malloc`
  and must be released by the caller via `cib_free_string`.
* Models live behind integer **handles** in a process-global registry; the
  caller passes the handle back on each call and releases it with `cib_free`.
* The surface is purely **integer-indexed**: scenarios are JSON arrays of
  0-based variant indices, and descriptor/variant references in
  `cib_set_impact` / `cib_get_impact` are 0-based indices. Human-readable
  name resolution is left to the caller (which receives the descriptor and
  variant names from `cib_load`).

Runtime lifecycle (from the compiled library): call the generated
`init_julia(argc, argv)` once before any `cib_*` call, and `shutdown_julia`
at the end. Pass `-t auto` (or `-tN`) in `argv` to enable threads.

Removed options
---------------
`cib_consistent` once accepted `seed`, `ignore_cycles` and `exhaustive`. All
three are gone, along with the simulated annealing and Monte-Carlo sampling
they configured. The engine now always enumerates the whole scenario space and
is fully deterministic, so there is no randomness to seed and no partial search
to opt out of. Callers passing those keys will simply have them ignored.
"""
module CIBCApi

using CrossImpactBalances

include("json.jl")   # parse_json_value, json_serialize, ... (dependency-free)

# ── Handle registry ─────────────────────────────────────────────────────────

const _MODELS = Dict{Cint,CIB}()
const _NEXT   = Ref{Cint}(Cint(1))
const _LOCK   = ReentrantLock()

function _register(cib::CIB)::Cint
    lock(_LOCK) do
        h = _NEXT[]
        _NEXT[] = h + Cint(1)
        _MODELS[h] = cib
        return h
    end
end

function _get(handle::Cint)::CIB
    lock(_LOCK) do
        haskey(_MODELS, handle) || error("Unknown model handle: $(Int(handle))")
        return _MODELS[handle]
    end
end

# ── String / JSON marshalling ───────────────────────────────────────────────

"""Copy a Julia `String` into malloc'd C memory and return it as a `Cstring`.
The caller frees it with `cib_free_string`."""
function _cstr(s::String)::Cstring
    n = sizeof(s)
    p = Libc.malloc(n + 1)
    p == C_NULL && throw(OutOfMemoryError())
    GC.@preserve s unsafe_copyto!(Ptr{UInt8}(p), pointer(s), n)
    unsafe_store!(Ptr{UInt8}(p), 0x00, n + 1)   # NUL terminator
    return Cstring(p)
end

_reply(obj)::Cstring = _cstr(json_serialize(obj))

function _ok(pairs::Pair...)::Cstring
    d = Dict{String,Any}("ok" => true)
    for (k, v) in pairs
        d[k] = v
    end
    return _reply(d)
end

_fail(msg::AbstractString)::Cstring =
    _reply(Dict{String,Any}("ok" => false, "error" => String(msg)))

"""Run `f()` and marshal its result; convert any exception into an error JSON
reply so exceptions never cross the C boundary."""
function _protect(f)::Cstring
    try
        return f()
    catch e
        return _fail(sprint(showerror, e))
    end
end

# Parse a JSON string arg into a Julia value.
_parse(arg::Cstring) = parse_json_value(unsafe_string(arg), 1)[1]

# A parsed JSON array of integers -> Vector{Int}.
_intvec(x)::Vector{Int} = Int[Int(v) for v in x]

# Optional field with a default from a parsed JSON object (or empty).
function _opt(obj, key::String, default)
    (obj isa AbstractDict && haskey(obj, key) && obj[key] !== nothing) ?
        obj[key] : default
end

_rule(name) = name == "sequential" ? SequentialSuccession() : GlobalSuccession()

# ── Exported C API ──────────────────────────────────────────────────────────

Base.@ccallable function cib_load(path::Cstring)::Cstring
    _protect() do
        p = unsafe_string(path)
        isfile(p) || return _fail("File not found: $p")
        cib = load_scw(p; kernel = Vector{Vector{Int}}())
        h = _register(cib)
        variants = Dict{String,Any}(d => cib.variants[d] for d in cib.descriptors)
        _ok("handle" => Int(h),
            "descriptors" => cib.descriptors,
            "variants" => variants,
            "n_descriptors" => cib.numberOfDescriptors,
            "n_scenarios" => max_signature(cib) + 1)
    end
end

Base.@ccallable function cib_free(handle::Cint)::Cvoid
    lock(_LOCK) do
        delete!(_MODELS, handle)
    end
    return nothing
end

Base.@ccallable function cib_copy(handle::Cint)::Cstring
    _protect() do
        cib = _get(handle)
        h = _register(deepcopy(cib))
        _ok("handle" => Int(h))
    end
end

Base.@ccallable function cib_consistent(handle::Cint, opts::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(opts)
        # Options are `{rule?, algorithm?}`. The search is always exhaustive and
        # always deterministic, so there is nothing else to configure — see the
        # note on the removed options in this module's docstring.
        kern = find_consistent(cib;
            rule = _rule(_opt(o, "rule", "global")),
            algorithm = Symbol(_opt(o, "algorithm", "auto")))
        _ok("scenarios" => kern)
    end
end

Base.@ccallable function cib_matrix(handle::Cint)::Cstring
    _protect() do
        cib = _get(handle)
        rows = [cib.cim[i, :] for i in 1:size(cib.cim, 1)]
        _ok("matrix" => rows)
    end
end

Base.@ccallable function cib_basins(handle::Cint, opts::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(opts)
        fps, sizes, cycles = find_basins(cib; rule = _rule(_opt(o, "rule", "global")))
        _ok("fixed_points" => fps,
            "basin_sizes" => sizes,
            "cycle_count" => cycles,
            "total" => max_signature(cib) + 1)
    end
end

Base.@ccallable function cib_impact_balance(handle::Cint, scenario::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        u = _intvec(_parse(scenario))
        _ok("ib" => impact_balance(cib, u))
    end
end

Base.@ccallable function cib_signature(handle::Cint, scenario::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        u = _intvec(_parse(scenario))
        _ok("signature" => signature(cib, u))
    end
end

Base.@ccallable function cib_inv_signature(handle::Cint, sig::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        s = Int(_parse(sig))
        _ok("scenario" => inv_signature(cib, s))
    end
end

Base.@ccallable function cib_succession(handle::Cint, req::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(req)
        u = _intvec(o["scenario"])
        rule = _rule(_opt(o, "rule", "global"))
        maxsteps = Int(_opt(o, "max_steps", max_signature(cib) + 10))

        steps = Vector{Vector{Int}}([copy(u)])
        seen = Set{Int}([signature(cib, u)])
        cycle_length = 0
        cur = copy(u)
        for _ in 1:maxsteps
            nxt = succession_step(rule, cib, cur)
            nsig = signature(cib, nxt)
            if nsig in seen
                if nxt == steps[end]
                    cycle_length = 1
                else
                    push!(steps, copy(nxt))
                    for k in length(steps):-1:1
                        cycle_length += 1
                        signature(cib, steps[k]) == nsig && break
                    end
                end
                break
            end
            push!(seen, nsig)
            push!(steps, copy(nxt))
            cur = nxt
        end
        _ok("steps" => steps, "cycle_length" => cycle_length)
    end
end

Base.@ccallable function cib_set_impact(handle::Cint, req::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(req)
        s = _intvec(o["source"])   # [descriptor_index, variant_index], 0-based
        t = _intvec(o["target"])
        v = Int(o["value"])
        old = set_impact!(cib, s[1], s[2], t[1], t[2], v)
        _ok("old" => old)
    end
end

Base.@ccallable function cib_get_impact(handle::Cint, req::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(req)
        s = _intvec(o["source"])
        t = _intvec(o["target"])
        _ok("value" => get_impact(cib, s[1], s[2], t[1], t[2]))
    end
end

Base.@ccallable function cib_free_string(ptr::Cstring)::Cvoid
    ptr == Cstring(C_NULL) || Libc.free(Ptr{Cvoid}(ptr))
    return nothing
end

end # module

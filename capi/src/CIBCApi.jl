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

Very large models
-----------------
A real ScenarioWizard matrix can have more scenarios than `Int64` can count —
10²⁴ is an ordinary size in the Weimer-Jehle corpus — and `find_consistent`
searches those anyway, because branch-and-bound never numbers a scenario. This
surface reaches them too: `cib_load` reports the space with `scenario_count`
(exact to ~1.7×10³⁸) rather than `max_signature + 1`, which refuses there, and
`cib_signature` answers with the `Int128` signature rather than the `Int64` one,
which wraps into negative numbers. Counts and signatures past 2⁵³ cross the
boundary as **decimal strings**, because a consumer that parses JSON numbers
into doubles — every JavaScript one, and anything using a float-backed parser —
silently loses their last digits. Python's `int()` accepts either form
unchanged; see `_json_count` below.

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

# ── Numbers too big for a JSON number ───────────────────────────────────────

"""Encode a count or a signature for JSON: a number while a `Float64` holds it
exactly, a decimal **string** past 2⁵³ so a consumer that parses numbers into
doubles cannot mangle its last digits. `int(...)` in Python takes either form
transparently. Mirrors `json_count` in `app/src/CIBApp.jl`."""
_json_count(n::Integer) = Int128(n) <= Int128(2)^53 ? Int(n) : string(n)

"""The inverse, for a signature arriving from a caller: accept the decimal
string `_json_count` emits as well as a plain number, and refuse anything
`inv_signature` cannot represent by name, rather than letting it surface as an
overflow from the JSON number parser."""
function _signature_arg(x)::Int
    signatureValue = x isa AbstractString ? parse(Int128, x) : Int128(x)
    signatureValue <= Int128(typemax(Int)) || throw(ArgumentError(
        "cib_inv_signature: signature $(signatureValue) is past typemax(Int64) = " *
        "$(typemax(Int)), and inv_signature works in Int64 — it cannot invert it. " *
        "(cib_signature still reports such signatures exactly; only the inverse " *
        "direction stops here.) The scenario's variant indices identify it at any " *
        "model size."))
    return Int(signatureValue)
end

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
            # scenario_count, not max_signature + 1: max_signature THROWS past
            # typemax(Int64), and this line is the whole reason a model that
            # size could not be loaded through the C API at all — the failure
            # came back from cib_load, so cib_consistent, which handles those
            # models perfectly well, was unreachable for them.
            "n_scenarios" => _json_count(scenario_count(cib)))
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
            "total" => _json_count(scenario_count(cib)))
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
        # _signature128, not signature: the same mixed-radix number, in a type
        # wide enough for every model cib_load now accepts. `signature` returns
        # an Int and wraps silently into negative numbers past typemax(Int64) —
        # unreachable while cib_load refused those models, and a wrong answer
        # rather than an error the moment it stopped refusing.
        _ok("signature" => _json_count(CrossImpactBalances._signature128(cib, u)))
    end
end

Base.@ccallable function cib_inv_signature(handle::Cint, sig::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        s = _signature_arg(_parse(sig))
        _ok("scenario" => inv_signature(cib, s))
    end
end

Base.@ccallable function cib_succession(handle::Cint, req::Cstring)::Cstring
    _protect() do
        cib = _get(handle)
        o = _parse(req)
        u = _intvec(o["scenario"])
        rule = _rule(_opt(o, "rule", "global"))
        # Int128 signatures as the "have I been here before" key, and a step
        # limit taken from scenario_count rather than max_signature: succession
        # itself costs the same at any model size, and cib_load now reaches
        # here with models whose Int64 signatures wrap — where two different
        # scenarios can share a key and the walk reports a cycle that is not
        # there. The limit only has to outlast the longest possible trajectory,
        # which cannot revisit a scenario; it is clamped so it stays an Int.
        # Mirrors do_succession in mcp/julia_worker.jl.
        maxsteps = Int(_opt(o, "max_steps",
                            Int(min(scenario_count(cib) + 10, Int128(typemax(Int))))))

        steps = Vector{Vector{Int}}([copy(u)])
        seen = Set{Int128}([CrossImpactBalances._signature128(cib, u)])
        cycle_length = 0
        cur = copy(u)
        for _ in 1:maxsteps
            nxt = succession_step(rule, cib, cur)
            nsig = CrossImpactBalances._signature128(cib, nxt)
            if nsig in seen
                if nxt == steps[end]
                    cycle_length = 1
                else
                    push!(steps, copy(nxt))
                    # Count back from the step BEFORE the one just pushed.
                    # Starting AT it matched immediately — `nxt` is what nsig
                    # is the signature of — so every cycle came back as
                    # length 1, which `Model.succession` in python/ turns into
                    # `converged: true`. CIB_global from [0,0,0] is a two-cycle
                    # that reported as a fixed point.
                    for k in (length(steps) - 1):-1:1
                        cycle_length += 1
                        CrossImpactBalances._signature128(cib, steps[k]) == nsig && break
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

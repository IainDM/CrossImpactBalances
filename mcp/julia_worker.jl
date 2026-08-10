"""
Persistent Julia worker process for the Node.js MCP server.

Loads CrossImpactBalances once, then sits in a loop reading commands
from stdin and writing JSON results to stdout.  This avoids paying
Julia's startup + compilation cost on every tool call.

Protocol (newline-delimited, tab-separated fields):
  Ready signal:  {"status":"ready"}
  Request:       op\tfile_path[\tDescriptor=Variant ...]
  Response:      {"ok":true,"text":"..."} or {"ok":false,"error":"..."}

Operations:
  fixed_points   — find all consistent scenarios (fixed points)
  basins         — exhaustive basin-of-attraction analysis
  succession     — trace succession from a starting scenario
  impact_balance — compute impact scores for a scenario

Not intended to be run directly — spawned by server.js.
"""

# ── Load package (one-time cost) ───────────────────────────────────────────

println(stderr, "julia_worker: loading CrossImpactBalances...")

using CrossImpactBalances

println(stderr, "julia_worker: package loaded")

# ── Single-entry CIB cache ────────────────────────────────────────────────

let
global _cached_path = ""
global _cached_cib = nothing
end

function get_cib(file_path::String)
    global _cached_path, _cached_cib
    if file_path != _cached_path
        _cached_cib = load_scw(file_path; kernel=Vector{Vector{Int}}())
        _cached_path = file_path
        println(stderr, "julia_worker: loaded $(basename(file_path))")
    end
    return _cached_cib
end

# ── Minimal JSON string escaping (output only) ────────────────────────────

function json_escape(s::AbstractString)
    buf = IOBuffer()
    for c in s
        if c == '"';       write(buf, "\\\"")
        elseif c == '\\';  write(buf, "\\\\")
        elseif c == '\n';  write(buf, "\\n")
        elseif c == '\r';  write(buf, "\\r")
        elseif c == '\t';  write(buf, "\\t")
        elseif codepoint(c) < 0x20
            write(buf, "\\u$(lpad(string(codepoint(c), base=16), 4, '0'))")
        else
            write(buf, c)
        end
    end
    return String(take!(buf))
end

function send_ok(text::String)
    println(stdout, "{\"ok\":true,\"text\":\"$(json_escape(text))\"}")
    flush(stdout)
    Base.Libc.flush_cstdio()
end

function send_err(msg::String)
    println(stdout, "{\"ok\":false,\"error\":\"$(json_escape(msg))\"}")
    flush(stdout)
    Base.Libc.flush_cstdio()
end

# ── Scenario name resolution ─────────────────────────────────────────────

"""
Parse descriptor=variant pairs and resolve to 0-based variant indices.
Returns a Vector{Int} of length ndesc, or throws on invalid names.
"""
function resolve_scenario(cib, pairs::Vector{SubString{String}})
    scenario = fill(-1, cib.numberOfDescriptors)
    for pair in pairs
        idx = findfirst('=', pair)
        isnothing(idx) && error("Invalid pair (missing '='): $pair")
        desc_name = String(pair[1:idx-1])
        var_name = String(pair[idx+1:end])

        desc_idx = findfirst(==(desc_name), cib.descriptors)
        isnothing(desc_idx) && error("Unknown descriptor: \"$desc_name\". Available: $(join(cib.descriptors, ", "))")

        vars = cib.variants[desc_name]
        var_idx = findfirst(==(var_name), vars)
        isnothing(var_idx) && error("Unknown variant \"$var_name\" for descriptor \"$desc_name\". Available: $(join(vars, ", "))")

        scenario[desc_idx] = var_idx - 1  # 0-based
    end

    for (i, desc) in enumerate(cib.descriptors)
        if scenario[i] == -1
            error("Missing descriptor: \"$desc\". All descriptors must be specified.")
        end
    end

    return scenario
end

# ── Result formatting ─────────────────────────────────────────────────────

function format_scenario(cib, u)
    lines = String[]
    for (j, desc) in enumerate(cib.descriptors)
        variant_name = cib.variants[desc][u[j] + 1]
        push!(lines, "  $desc = $variant_name")
    end
    return join(lines, "\n")
end

function format_model_header(buf, cib)
    println(buf, "Descriptors: $(cib.numberOfDescriptors)")
    print(buf, "Variants per descriptor: [")
    print(buf, join(cib.numberOfVariants, ", "))
    println(buf, "]")
    # scenario_count, not max_signature + 1: this header is printed alongside
    # results that find_consistent can now produce for models whose signatures
    # overrun Int64, and max_signature would throw here and discard them.
    total = scenario_count(cib)
    println(buf, "Total scenario space: $total")
    println(buf)

    println(buf, "Descriptors and variants:")
    for desc in cib.descriptors
        vars = cib.variants[desc]
        println(buf, "  $desc: $(join(vars, ", "))")
    end
    println(buf)
end

function format_fixed_points(cib, fps)
    buf = IOBuffer()

    println(buf, "Consistent Scenarios (Fixed Points)")
    println(buf, "=" ^ 50)
    println(buf)

    format_model_header(buf, cib)

    println(buf, "Fixed Points Found: $(length(fps))")
    println(buf, "-" ^ 50)

    for (rank, fp) in enumerate(fps)
        println(buf)
        println(buf, "Scenario #$rank")
        println(buf, format_scenario(cib, fp))
    end

    return String(take!(buf))
end

function format_basins(cib, fps, basins, cycle_count, total)
    buf = IOBuffer()

    println(buf, "Basin-of-Attraction Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    format_model_header(buf, cib)

    perm = sortperm(basins, rev=true)
    println(buf, "Consistent Scenarios (Fixed Points): $(length(fps))")
    println(buf, "-" ^ 50)

    for (rank, idx) in enumerate(perm)
        fp = fps[idx]
        basin = basins[idx]
        pct = round(100.0 * basin / total; digits=2)

        println(buf)
        println(buf, "Scenario #$rank (basin: $basin scenarios, $pct%)")
        println(buf, format_scenario(cib, fp))
    end

    if cycle_count > 0
        println(buf)
        cyc_pct = round(100.0 * cycle_count / total; digits=2)
        println(buf, "Scenarios in cycles (non-convergent): $cycle_count ($cyc_pct%)")
    end

    println(buf)
    converging = sum(basins)
    println(buf, "Coverage: $converging converging + $cycle_count cycling = $total total")

    return String(take!(buf))
end

function format_succession(cib, steps, cycle_length)
    buf = IOBuffer()

    println(buf, "Succession Trace")
    println(buf, "=" ^ 50)
    println(buf)

    for (i, step) in enumerate(steps)
        if i == 1
            println(buf, "Start:")
        else
            println(buf, "Step $(i-1):")
        end
        println(buf, format_scenario(cib, step))
        println(buf)
    end

    if cycle_length == 1
        println(buf, "Result: converges to a fixed point (consistent scenario) in $(length(steps) - 1) step(s)")
    else
        println(buf, "Result: enters a cycle of length $cycle_length (no fixed point)")
    end

    return String(take!(buf))
end

function format_impact_balance(cib, scenario, ib)
    buf = IOBuffer()

    println(buf, "Impact Balance Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    println(buf, "Scenario:")
    println(buf, format_scenario(cib, scenario))
    println(buf)

    println(buf, "Impact scores per variant:")
    idx = 1
    is_fixed = true
    for (j, desc) in enumerate(cib.descriptors)
        nv = cib.numberOfVariants[j]
        scores = ib[idx:idx + nv - 1]
        max_score = maximum(scores)
        selected_score = scores[scenario[j] + 1]
        if selected_score < max_score
            is_fixed = false
        end

        println(buf, "  $desc:")
        for (k, var) in enumerate(cib.variants[desc])
            score = scores[k]
            marker = ""
            if k - 1 == scenario[j]
                marker = " [selected]"
            end
            if score == max_score && max_score > selected_score && k - 1 != scenario[j]
                marker *= " ← would win"
            end
            println(buf, "    $var: $score$marker")
        end
        idx += nv
    end

    println(buf)
    if is_fixed
        println(buf, "This scenario IS self-consistent (fixed point).")
        println(buf, "All selected variants have the highest impact scores in their descriptor.")
    else
        println(buf, "This scenario is NOT self-consistent.")
        println(buf, "Under succession, some variants would change (see \"would win\" above).")
    end

    return String(take!(buf))
end

# ── Operation handlers ────────────────────────────────────────────────────

function do_fixed_points(file_path::String)
    cib = get_cib(file_path)
    fps = find_consistent(cib)
    text = format_fixed_points(cib, fps)
    println(stderr, "julia_worker: done ($(length(fps)) fixed points)")
    send_ok(text)
end

function do_basins(file_path::String)
    cib = get_cib(file_path)
    fps, basins, cycle_count = find_basins(cib)
    total = max_signature(cib) + 1
    text = format_basins(cib, fps, basins, cycle_count, total)
    println(stderr, "julia_worker: done ($(length(fps)) fixed points, basins computed)")
    send_ok(text)
end

function do_succession(file_path::String, pairs::Vector{SubString{String}})
    cib = get_cib(file_path)
    scenario = resolve_scenario(cib, pairs)

    # Collect steps. Scenarios are keyed by their Int128 signature, not the
    # Int64 one: succession works at any model size, and load_scw now reaches
    # here with models whose Int64 signatures wrap — where two different
    # scenarios can share a key and the walk would report a cycle that is not
    # there. The step limit is the scenario count for the same reason (a
    # trajectory cannot outlast the space), capped so it stays an Int.
    steps = [copy(scenario)]
    v = copy(scenario)
    seen_sigs = Set{Int128}([CrossImpactBalances._signature128(cib, v)])
    cycle_length = 0
    step_limit = Int(min(scenario_count(cib) + 10, Int128(typemax(Int))))

    for _ in 1:step_limit
        v = succession_step(cib, v)
        v_sig = CrossImpactBalances._signature128(cib, v)
        if v_sig in seen_sigs
            # Check if it's a fixed point (same as last step) or a cycle
            if v == steps[end]
                cycle_length = 1
            else
                push!(steps, copy(v))
                # Count cycle length
                cycle_length = 0
                for k in length(steps):-1:1
                    cycle_length += 1
                    if CrossImpactBalances._signature128(cib, steps[k]) == v_sig
                        break
                    end
                end
            end
            break
        end
        push!(seen_sigs, v_sig)
        push!(steps, copy(v))
    end

    text = format_succession(cib, steps, cycle_length)
    println(stderr, "julia_worker: succession traced ($(length(steps)) steps, cycle=$(cycle_length))")
    send_ok(text)
end

function do_impact_balance(file_path::String, pairs::Vector{SubString{String}})
    cib = get_cib(file_path)
    scenario = resolve_scenario(cib, pairs)

    ib = impact_balance(cib, scenario)
    text = format_impact_balance(cib, scenario, ib)
    println(stderr, "julia_worker: impact balance computed")
    send_ok(text)
end

# ── Main loop ─────────────────────────────────────────────────────────────

function main()
    # Signal readiness to Node.js
    println(stdout, "{\"status\":\"ready\"}")
    flush(stdout)
    Base.Libc.flush_cstdio()

    println(stderr, "julia_worker: ready, waiting for commands")

    for line in eachline(stdin)
        stripped = strip(line)
        isempty(stripped) && continue

        parts = split(stripped, '\t')
        op = parts[1]

        println(stderr, "julia_worker: op=$op")

        try
            if length(parts) < 2
                send_err("Missing file_path argument")
                continue
            end
            file_path = String(parts[2])

            if !isfile(file_path)
                send_err("File not found: $file_path")
                continue
            end

            if op == "fixed_points"
                do_fixed_points(file_path)
            elseif op == "basins"
                do_basins(file_path)
            elseif op == "succession"
                if length(parts) < 3
                    send_err("Missing scenario argument (Descriptor=Variant pairs)")
                    continue
                end
                do_succession(file_path, parts[3:end])
            elseif op == "impact_balance"
                if length(parts) < 3
                    send_err("Missing scenario argument (Descriptor=Variant pairs)")
                    continue
                end
                do_impact_balance(file_path, parts[3:end])
            else
                send_err("Unknown operation: $op")
            end
        catch e
            msg = "Error: $(sprint(showerror, e))"
            println(stderr, "julia_worker: $msg")
            send_err(msg)
        end
    end

    println(stderr, "julia_worker: shutting down (EOF)")
end

main()

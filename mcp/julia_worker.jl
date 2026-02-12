"""
Persistent Julia worker process for the Node.js MCP server.

Loads CrossImpactBalances once, then sits in a loop reading file paths
from stdin and writing JSON results to stdout.  This avoids paying
Julia's startup + compilation cost on every tool call.

Protocol (newline-delimited):
  Ready signal:  {"status":"ready"}
  Request:       <file_path>          (one line)
  Response:      {"ok":true,"text":"..."} or {"ok":false,"error":"..."}

Not intended to be run directly — spawned by server.js.
"""

# ── Load package (one-time cost) ───────────────────────────────────────────

println(stderr, "julia_worker: loading CrossImpactBalances...")

using CrossImpactBalances

println(stderr, "julia_worker: package loaded")

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

# ── Result formatting ─────────────────────────────────────────────────────

function format_results(cib, fps, basins, cycle_count, total)
    buf = IOBuffer()

    println(buf, "Cross-Impact Balance Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    println(buf, "Descriptors: $(cib.ndesc)")
    print(buf, "Variants per descriptor: [")
    print(buf, join(cib.nvariants, ", "))
    println(buf, "]")
    println(buf, "Total scenario space: $total")
    println(buf)

    println(buf, "Descriptors and variants:")
    for desc in cib.descriptors
        vars = cib.variants[desc]
        println(buf, "  $desc: $(join(vars, ", "))")
    end
    println(buf)

    perm = sortperm(basins, rev=true)
    println(buf, "Consistent Scenarios (Fixed Points): $(length(fps))")
    println(buf, "-" ^ 50)

    for (rank, idx) in enumerate(perm)
        fp = fps[idx]
        basin = basins[idx]
        pct = round(100.0 * basin / total; digits=2)

        println(buf)
        println(buf, "Scenario #$rank (basin: $basin scenarios, $pct%)")

        for (j, desc) in enumerate(cib.descriptors)
            variant_name = cib.variants[desc][fp[j] + 1]
            println(buf, "  $desc = $variant_name")
        end
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

# ── Main loop ─────────────────────────────────────────────────────────────

function main()
    # Signal readiness to Node.js
    println(stdout, "{\"status\":\"ready\"}")
    flush(stdout)
    Base.Libc.flush_cstdio()

    println(stderr, "julia_worker: ready, waiting for commands")

    for line in eachline(stdin)
        file_path = strip(line)
        isempty(file_path) && continue

        println(stderr, "julia_worker: analyzing $(basename(file_path))")

        try
            if !isfile(file_path)
                send_err("File not found: $file_path")
                continue
            end

            cib = load_scw(file_path; kernel=Vector{Vector{Int}}())
            fps, basins, cycle_count = find_basins(cib)
            total = max_signature(cib) + 1
            text = format_results(cib, fps, basins, cycle_count, total)

            println(stderr, "julia_worker: done ($(length(fps)) fixed points)")
            send_ok(text)
        catch e
            msg = "Analysis failed: $(sprint(showerror, e))"
            println(stderr, "julia_worker: error: $msg")
            send_err(msg)
        end
    end

    println(stderr, "julia_worker: shutting down (EOF)")
end

main()

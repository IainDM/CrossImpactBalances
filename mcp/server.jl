"""
MCP (Model Context Protocol) server for CrossImpactBalances.

Exposes a `scw_fixed_points` tool that accepts a ScenarioWizard .scw file
and returns the fixed points (consistent scenarios) together with their
basin-of-attraction sizes.

Run directly:
    julia --project=/path/to/CrossImpactBalances /path/to/CrossImpactBalances/mcp/server.jl

Add to Claude Code (recommended — uses run-server.sh wrapper so that Julia
is found even when spawned by a non-interactive host process):
    claude mcp add crossimpactbalances -- /path/to/CrossImpactBalances/mcp/run-server.sh
"""

# Lazy-load CrossImpactBalances on first tool call, not at startup.
# The MCP handshake (initialize, tools/list) needs no package symbols and
# must complete quickly.  Package compilation is deferred to the first
# tools/call, whose timeout is much longer.
const _cib_loaded = Ref(false)

function _ensure_cib_loaded()
    _cib_loaded[] && return
    @eval begin
        using CrossImpactBalances
        function _run_analysis(file_path::AbstractString)
            cib = load_scw(file_path; kernel=Vector{Vector{Int}}())
            fps, basins, cycle_count = find_basins(cib)
            total = max_signature(cib) + 1
            return (cib, fps, basins, cycle_count, total)
        end
    end
    _cib_loaded[] = true
end

# ── Minimal JSON parser/serializer (shared with test/benchmark.jl) ──────────

include(joinpath(@__DIR__, "json.jl"))

# ── MCP transport (Content-Length framing over stdio) ────────────────────────

function read_message(input::IO)
    content_length = 0
    while true
        line = readline(input)
        stripped = strip(line)
        if isempty(stripped)
            break
        end
        if startswith(stripped, "Content-Length:")
            content_length = parse(Int, strip(stripped[16:end]))
        end
    end
    if content_length == 0
        error("Missing Content-Length header")
    end
    body = read(input, content_length)
    return parse_json_value(String(body), 1)[1]
end

function send_message(output::IO, msg::Dict)
    body = json_serialize(msg)
    header = "Content-Length: $(sizeof(body))\r\n\r\n"
    write(output, header)
    write(output, body)
    flush(output)
    Base.Libc.flush_cstdio()  # also flush C-level stdout buffer (pipe buffering)
end

# ── MCP protocol handlers ───────────────────────────────────────────────────

const SERVER_NAME = "crossimpactbalances-mcp"
const SERVER_VERSION = "0.1.0"
const PROTOCOL_VERSION = "2024-11-05"

function handle_initialize(id)
    return Dict(
        "jsonrpc" => "2.0",
        "id" => id,
        "result" => Dict(
            "protocolVersion" => PROTOCOL_VERSION,
            "capabilities" => Dict(
                "tools" => Dict()
            ),
            "serverInfo" => Dict(
                "name" => SERVER_NAME,
                "version" => SERVER_VERSION
            )
        )
    )
end

function handle_tools_list(id)
    return Dict(
        "jsonrpc" => "2.0",
        "id" => id,
        "result" => Dict(
            "tools" => [
                Dict(
                    "name" => "scw_fixed_points",
                    "description" => "Analyze a ScenarioWizard (.scw) cross-impact balance file. " *
                        "Finds all consistent scenarios (fixed points of the succession " *
                        "operator) and their basins of attraction: the number of starting " *
                        "combinations that converge to each fixed point.",
                    "inputSchema" => Dict(
                        "type" => "object",
                        "properties" => Dict(
                            "file_path" => Dict(
                                "type" => "string",
                                "description" => "Absolute path to a ScenarioWizard .scw file"
                            ),
                            "scw_content" => Dict(
                                "type" => "string",
                                "description" => "Raw contents of a ScenarioWizard .scw file. " *
                                    "Use this when the file content is available directly " *
                                    "(e.g. uploaded or pasted) instead of as a filesystem path. " *
                                    "Provide either file_path or scw_content, not both."
                            )
                        )
                    )
                )
            ]
        )
    )
end

function handle_tools_call(id, params)
    name = get(params, "name", "")
    args = get(params, "arguments", Dict())

    if name == "scw_fixed_points"
        return execute_scw_fixed_points(id, args)
    else
        return make_error(id, -32601, "Unknown tool: $name")
    end
end

# ── Tool implementation ─────────────────────────────────────────────────────

function execute_scw_fixed_points(id, args)
    file_path = get(args, "file_path", "")
    scw_content = get(args, "scw_content", "")

    if isempty(file_path) && isempty(scw_content)
        return make_error(id, -32602, "Either file_path or scw_content must be provided")
    end

    if !isempty(file_path) && !isempty(scw_content)
        return make_error(id, -32602, "Provide either file_path or scw_content, not both")
    end

    temp_file = ""
    if !isempty(scw_content)
        temp_file = tempname() * ".scw"
        write(temp_file, scw_content)
        file_path = temp_file
    end

    if !isfile(file_path)
        !isempty(temp_file) && isfile(temp_file) && rm(temp_file)
        return make_error(id, -32602, "File not found: $file_path")
    end

    try
        _ensure_cib_loaded()
        cib, fps, basins, cycle_count, total = Base.invokelatest(_run_analysis, file_path)
        text = format_results(cib, fps, basins, cycle_count, total)

        return Dict(
            "jsonrpc" => "2.0",
            "id" => id,
            "result" => Dict(
                "content" => [
                    Dict("type" => "text", "text" => text)
                ]
            )
        )
    catch e
        return make_error(id, -32603, "Analysis failed: $(sprint(showerror, e))")
    finally
        !isempty(temp_file) && isfile(temp_file) && rm(temp_file)
    end
end

function format_results(cib, fps, basins, cycle_count, total)
    buf = IOBuffer()

    println(buf, "Cross-Impact Balance Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    # Problem size
    println(buf, "Descriptors: $(cib.numberOfDescriptors)")
    print(buf, "Variants per descriptor: [")
    print(buf, join(cib.numberOfVariants, ", "))
    println(buf, "]")
    println(buf, "Total scenario space: $total")
    println(buf)

    # Descriptor listing
    println(buf, "Descriptors and variants:")
    for desc in cib.descriptors
        vars = cib.variants[desc]
        println(buf, "  $desc: $(join(vars, ", "))")
    end
    println(buf)

    # Fixed points sorted by basin size (descending)
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
            variant_name = cib.variants[desc][fp[j] + 1]  # 0-based to 1-based
            println(buf, "  $desc = $variant_name")
        end
    end

    # Cycles
    if cycle_count > 0
        println(buf)
        cyc_pct = round(100.0 * cycle_count / total; digits=2)
        println(buf, "Scenarios in cycles (non-convergent): $cycle_count ($cyc_pct%)")
    end

    # Summary
    println(buf)
    converging = sum(basins)
    println(buf, "Coverage: $converging converging + $cycle_count cycling = $total total")

    return String(take!(buf))
end

# ── Helpers ──────────────────────────────────────────────────────────────────

function make_error(id, code, message)
    return Dict(
        "jsonrpc" => "2.0",
        "id" => id,
        "error" => Dict(
            "code" => code,
            "message" => message
        )
    )
end

# ── Main event loop ─────────────────────────────────────────────────────────

function main()
    println(stderr, "crossimpactbalances-mcp: starting (Julia ", VERSION, ")")

    # Try to disable C-level stdout buffering (helps on some platforms).
    # Wrapped in try-catch because the C symbol lookup can fail on some
    # Windows configurations; the flush() calls in send_message are the
    # primary mechanism and always work.
    try
        if Sys.iswindows()
            cstdout_ptr = ccall((:__acrt_iob_func, "ucrtbase"), Ptr{Cvoid}, (UInt32,), 1)
            ionbf = 0x0004  # _IONBF on Windows UCRT
        else
            cstdout_ptr = unsafe_load(cglobal(:stdout, Ptr{Cvoid}))
            ionbf = 2  # _IONBF on POSIX
        end
        ccall(:setvbuf, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Csize_t),
              cstdout_ptr, C_NULL, ionbf, 0)
    catch e
        println(stderr, "crossimpactbalances-mcp: setvbuf failed (non-fatal): ",
                sprint(showerror, e))
    end

    input = stdin
    output = stdout

    println(stderr, "crossimpactbalances-mcp: ready, waiting for messages")

    while !eof(input)
        local msg
        try
            msg = read_message(input)
        catch e
            if isa(e, EOFError)
                break
            end
            println(stderr, "crossimpactbalances-mcp: read error: ",
                    sprint(showerror, e))
            continue
        end

        method = get(msg, "method", nothing)
        id = get(msg, "id", nothing)
        params = get(msg, "params", Dict())

        println(stderr, "crossimpactbalances-mcp: ← ", something(method, "?"),
                !isnothing(id) ? " (id=$id)" : "")

        response = nothing

        if method == "initialize"
            response = handle_initialize(id)
        elseif method == "notifications/initialized"
            # Notification — no response
        elseif method == "tools/list"
            response = handle_tools_list(id)
        elseif method == "tools/call"
            response = handle_tools_call(id, params)
        elseif !isnothing(id)
            response = make_error(id, -32601, "Method not found: $method")
        end

        if !isnothing(response)
            send_message(output, response)
            println(stderr, "crossimpactbalances-mcp: → response for ",
                    something(method, "?"), " (id=$id)")
        end
    end

    println(stderr, "crossimpactbalances-mcp: shutting down (EOF)")
end

main()

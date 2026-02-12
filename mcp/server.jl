"""
MCP (Model Context Protocol) server for CrossImpactBalances.

Exposes a `scw_fixed_points` tool that accepts a ScenarioWizard .scw file
and returns the fixed points (consistent scenarios) together with their
basin-of-attraction sizes.

Run:
    julia --project=/path/to/CrossImpactBalances /path/to/CrossImpactBalances/mcp/server.jl

Add to Claude Code:
    claude mcp add crossimpactbalances -- julia --project=/path/to/CrossImpactBalances \\
        /path/to/CrossImpactBalances/mcp/server.jl
"""

using CrossImpactBalances

# ── Minimal JSON parser (from test/benchmark.jl) ────────────────────────────

function skip_ws(s, i)
    while i <= length(s) && s[i] in (' ', '\t', '\n', '\r')
        i += 1
    end
    return i
end

function parse_json_value(s, i)
    i = skip_ws(s, i)
    i > length(s) && error("Unexpected end of JSON")
    c = s[i]
    if c == '"'
        return parse_json_string(s, i)
    elseif c == '{'
        return parse_json_object(s, i)
    elseif c == '['
        return parse_json_array(s, i)
    elseif c == 't'
        return (true, i + 4)
    elseif c == 'f'
        return (false, i + 5)
    elseif c == 'n'
        return (nothing, i + 4)
    elseif isdigit(c) || c == '-'
        return parse_json_number(s, i)
    else
        error("Unexpected character at position $i: '$(s[i])'")
    end
end

function parse_json_string(s, i)
    i += 1  # skip opening "
    buf = IOBuffer()
    while i <= length(s) && s[i] != '"'
        if s[i] == '\\'
            i += 1
            if s[i] == 'n'; write(buf, '\n')
            elseif s[i] == 't'; write(buf, '\t')
            elseif s[i] == 'r'; write(buf, '\r')
            elseif s[i] == '"'; write(buf, '"')
            elseif s[i] == '\\'; write(buf, '\\')
            elseif s[i] == '/'; write(buf, '/')
            elseif s[i] == 'u'
                # \uXXXX Unicode escape
                hex = s[i+1:i+4]
                cp = parse(UInt32, hex; base=16)
                i += 4
                # Handle surrogate pairs
                if 0xD800 <= cp <= 0xDBFF && i + 2 <= length(s) && s[i+1] == '\\' && s[i+2] == 'u'
                    hex2 = s[i+3:i+6]
                    lo = parse(UInt32, hex2; base=16)
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00)
                    i += 6
                end
                write(buf, Char(cp))
            else write(buf, s[i])
            end
        else
            write(buf, s[i])
        end
        i += 1
    end
    return (String(take!(buf)), i + 1)  # skip closing "
end

function parse_json_number(s, i)
    j = i
    if j <= length(s) && s[j] == '-'
        j += 1
    end
    while j <= length(s) && isdigit(s[j])
        j += 1
    end
    is_float = false
    if j <= length(s) && s[j] == '.'
        is_float = true
        j += 1
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    if j <= length(s) && s[j] in ('e', 'E')
        is_float = true
        j += 1
        if j <= length(s) && s[j] in ('+', '-')
            j += 1
        end
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    numstr = s[i:j-1]
    val = is_float ? parse(Float64, numstr) : parse(Int, numstr)
    return (val, j)
end

function parse_json_array(s, i)
    i += 1  # skip [
    arr = Any[]
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == ']'
        return (arr, i + 1)
    end
    while true
        val, i = parse_json_value(s, i)
        push!(arr, val)
        i = skip_ws(s, i)
        if i > length(s) || s[i] == ']'
            return (arr, i + 1)
        end
        i += 1  # skip comma
    end
end

function parse_json_object(s, i)
    i += 1  # skip {
    obj = Dict{String, Any}()
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == '}'
        return (obj, i + 1)
    end
    while true
        i = skip_ws(s, i)
        key, i = parse_json_string(s, i)
        i = skip_ws(s, i)
        i += 1  # skip colon
        val, i = parse_json_value(s, i)
        obj[key] = val
        i = skip_ws(s, i)
        if i > length(s) || s[i] == '}'
            return (obj, i + 1)
        end
        i += 1  # skip comma
    end
end

# ── Minimal JSON serializer ─────────────────────────────────────────────────

function json_serialize(v::AbstractString)
    buf = IOBuffer()
    write(buf, '"')
    for c in v
        if c == '"';       write(buf, "\\\"")
        elseif c == '\\';  write(buf, "\\\\")
        elseif c == '\n';  write(buf, "\\n")
        elseif c == '\r';  write(buf, "\\r")
        elseif c == '\t';  write(buf, "\\t")
        elseif codepoint(c) < 0x20 || codepoint(c) >= 0x80
            # Control characters and non-ASCII → \uXXXX
            cp = codepoint(c)
            if cp <= 0xFFFF
                write(buf, "\\u$(lpad(string(cp, base=16), 4, '0'))")
            else
                # Surrogate pair for supplementary plane
                cp -= 0x10000
                hi = 0xD800 + (cp >> 10)
                lo = 0xDC00 + (cp & 0x3FF)
                write(buf, "\\u$(lpad(string(hi, base=16), 4, '0'))")
                write(buf, "\\u$(lpad(string(lo, base=16), 4, '0'))")
            end
        else
            write(buf, c)
        end
    end
    write(buf, '"')
    return String(take!(buf))
end

json_serialize(v::Integer) = string(v)
json_serialize(v::AbstractFloat) = string(v)
json_serialize(v::Bool) = v ? "true" : "false"
json_serialize(::Nothing) = "null"

function json_serialize(v::AbstractVector)
    return "[" * join([json_serialize(x) for x in v], ",") * "]"
end

function json_serialize(v::Tuple)
    return "[" * join([json_serialize(x) for x in v], ",") * "]"
end

function json_serialize(v::AbstractDict)
    pairs = [json_serialize(string(k)) * ":" * json_serialize(val) for (k, val) in v]
    return "{" * join(pairs, ",") * "}"
end

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
                            )
                        ),
                        "required" => ["file_path"]
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

    if isempty(file_path)
        return make_error(id, -32602, "Missing required parameter: file_path")
    end

    if !isfile(file_path)
        return make_error(id, -32602, "File not found: $file_path")
    end

    try
        # Load SCW with empty kernel (find_basins will compute everything)
        cib = load_scw(file_path; kernel=Vector{Vector{Int}}())
        fps, basins, cycle_count = find_basins(cib)
        text = format_results(cib, fps, basins, cycle_count)

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
    end
end

function format_results(cib, fps, basins, cycle_count)
    total = max_signature(cib) + 1
    buf = IOBuffer()

    println(buf, "Cross-Impact Balance Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    # Problem size
    println(buf, "Descriptors: $(cib.ndesc)")
    print(buf, "Variants per descriptor: [")
    print(buf, join(cib.nvariants, ", "))
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
    # Disable C-level stdout buffering — critical for MCP over pipes.
    # Julia buffers stdout when connected to a pipe (not a TTY),
    # so the MCP host may never see responses without this.
    let cstdout = cglobal(:stdout, Ptr{Cvoid})
        ccall(:setvbuf, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Csize_t),
              unsafe_load(cstdout), C_NULL, 2, 0)  # 2 = _IONBF (unbuffered)
    end

    input = stdin
    output = stdout

    while !eof(input)
        local msg
        try
            msg = read_message(input)
        catch e
            if isa(e, EOFError)
                break
            end
            println(stderr, "Error reading message: ", sprint(showerror, e))
            continue
        end

        method = get(msg, "method", nothing)
        id = get(msg, "id", nothing)
        params = get(msg, "params", Dict())

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
        end
    end
end

main()

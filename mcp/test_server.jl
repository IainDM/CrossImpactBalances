"""
Integration test for the MCP server.

Spawns the server as a subprocess and sends MCP protocol messages to validate
that initialization, tool listing, and tool execution work correctly.

Run:
    julia --project=. mcp/test_server.jl
"""

using Test

const PROJECT_DIR = dirname(@__DIR__)
const SERVER_SCRIPT = joinpath(@__DIR__, "server.jl")
const SAMPLE_DIR = joinpath(PROJECT_DIR, "test", "sample_files")

# ── Minimal JSON parser (duplicated here so tests are self-contained) ────────

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
    if c == '"';       return parse_json_string(s, i)
    elseif c == '{';   return parse_json_object(s, i)
    elseif c == '[';   return parse_json_array(s, i)
    elseif c == 't';   return (true, i + 4)
    elseif c == 'f';   return (false, i + 5)
    elseif c == 'n';   return (nothing, i + 4)
    elseif isdigit(c) || c == '-'; return parse_json_number(s, i)
    else error("Unexpected character at position $i: '$(s[i])'")
    end
end

function parse_json_string(s, i)
    i += 1
    buf = IOBuffer()
    while i <= length(s) && s[i] != '"'
        if s[i] == '\\'
            i += 1
            if s[i] == 'n'; write(buf, '\n')
            elseif s[i] == 't'; write(buf, '\t')
            elseif s[i] == 'r'; write(buf, '\r')
            elseif s[i] == '"'; write(buf, '"')
            elseif s[i] == '\\'; write(buf, '\\')
            elseif s[i] == 'u'
                hex = s[i+1:i+4]
                cp = parse(UInt32, hex; base=16)
                i += 4
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
    return (String(take!(buf)), i + 1)
end

function parse_json_number(s, i)
    j = i
    if j <= length(s) && s[j] == '-'; j += 1; end
    while j <= length(s) && isdigit(s[j]); j += 1; end
    is_float = false
    if j <= length(s) && s[j] == '.'
        is_float = true; j += 1
        while j <= length(s) && isdigit(s[j]); j += 1; end
    end
    if j <= length(s) && s[j] in ('e', 'E')
        is_float = true; j += 1
        if j <= length(s) && s[j] in ('+', '-'); j += 1; end
        while j <= length(s) && isdigit(s[j]); j += 1; end
    end
    val = is_float ? parse(Float64, s[i:j-1]) : parse(Int, s[i:j-1])
    return (val, j)
end

function parse_json_array(s, i)
    i += 1; arr = Any[]; i = skip_ws(s, i)
    if i <= length(s) && s[i] == ']'; return (arr, i + 1); end
    while true
        val, i = parse_json_value(s, i); push!(arr, val); i = skip_ws(s, i)
        if i > length(s) || s[i] == ']'; return (arr, i + 1); end
        i += 1
    end
end

function parse_json_object(s, i)
    i += 1; obj = Dict{String,Any}(); i = skip_ws(s, i)
    if i <= length(s) && s[i] == '}'; return (obj, i + 1); end
    while true
        i = skip_ws(s, i); key, i = parse_json_string(s, i)
        i = skip_ws(s, i); i += 1
        val, i = parse_json_value(s, i); obj[key] = val; i = skip_ws(s, i)
        if i > length(s) || s[i] == '}'; return (obj, i + 1); end
        i += 1
    end
end

# ── MCP message helpers ─────────────────────────────────────────────────────

function mcp_frame(json_str::String)
    return "Content-Length: $(sizeof(json_str))\r\n\r\n$json_str"
end

function read_mcp_response(io::IO)
    content_length = 0
    while true
        line = readline(io)
        stripped = strip(line)
        if isempty(stripped)
            break
        end
        if startswith(stripped, "Content-Length:")
            content_length = parse(Int, strip(stripped[16:end]))
        end
    end
    body = read(io, content_length)
    return parse_json_value(String(body), 1)[1]
end

# ── Tests ────────────────────────────────────────────────────────────────────

@testset "MCP Server" begin
    julia_cmd = Base.julia_cmd()
    cmd = `$julia_cmd --project=$PROJECT_DIR $SERVER_SCRIPT`
    proc = open(cmd, "r+")

    try
        # ── 1. Initialize ──
        @testset "initialize" begin
            init_msg = """{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}"""
            write(proc, mcp_frame(init_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 1
            @test haskey(resp, "result")
            @test resp["result"]["protocolVersion"] == "2024-11-05"
            @test haskey(resp["result"]["capabilities"], "tools")
            @test resp["result"]["serverInfo"]["name"] == "crossimpactbalances-mcp"
        end

        # ── 2. Initialized notification (no response expected) ──
        notif_msg = """{"jsonrpc":"2.0","method":"notifications/initialized"}"""
        write(proc, mcp_frame(notif_msg))
        flush(proc)

        # ── 3. tools/list ──
        @testset "tools/list" begin
            list_msg = """{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}"""
            write(proc, mcp_frame(list_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 2
            tools = resp["result"]["tools"]
            @test length(tools) == 1
            @test tools[1]["name"] == "scw_fixed_points"
            @test haskey(tools[1], "inputSchema")
            schema = tools[1]["inputSchema"]
            @test "file_path" in schema["required"]
        end

        # ── 4. tools/call with CIB_global.scw ──
        @testset "tools/call — CIB_global.scw" begin
            scw_path = joinpath(SAMPLE_DIR, "CIB_global.scw")
            # Escape backslashes in path for JSON (relevant on Windows)
            escaped_path = replace(scw_path, "\\" => "\\\\")
            call_msg = """{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"scw_fixed_points","arguments":{"file_path":"$escaped_path"}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 3
            @test haskey(resp, "result")
            content = resp["result"]["content"]
            @test length(content) >= 1
            @test content[1]["type"] == "text"

            text = content[1]["text"]

            # Must mention all 4 fixed points' variant names
            @test occursin("Consistent Scenarios (Fixed Points): 4", text)

            # Descriptor names must appear
            @test occursin("WTRD", text)
            @test occursin("WSEC", text)
            @test occursin("WECO", text)

            # Known variant assignments must appear
            @test occursin("Ntl", text)
            @test occursin("Mix", text)
            @test occursin("Alrt", text)
            @test occursin("ModGr", text)

            # Total scenario space
            @test occursin("Total scenario space: 36", text)

            # Coverage line must balance
            @test occursin("= 36 total", text)
        end

        # ── 5. tools/call with non-existent file ──
        @testset "tools/call — missing file" begin
            call_msg = """{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"scw_fixed_points","arguments":{"file_path":"/tmp/nonexistent.scw"}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 4
            @test haskey(resp, "error")
            @test occursin("not found", resp["error"]["message"])
        end

        # ── 6. Unknown tool ──
        @testset "tools/call — unknown tool" begin
            call_msg = """{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"bogus_tool","arguments":{}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 5
            @test haskey(resp, "error")
            @test occursin("Unknown tool", resp["error"]["message"])
        end

    finally
        close(proc)
    end
end

println("\nAll MCP server tests passed.")

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

# ── Minimal JSON parser (shared with mcp/server.jl and test/runtests.jl) ────

include(joinpath(@__DIR__, "json.jl"))

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
            @test haskey(schema["properties"], "file_path")
            @test haskey(schema["properties"], "scw_content")
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

        # ── 5. tools/call with scw_content ──
        @testset "tools/call — scw_content" begin
            scw_path = joinpath(SAMPLE_DIR, "CIB_global.scw")
            scw_content = read(scw_path, String)
            # Escape for JSON embedding
            escaped_content = replace(scw_content, "\\" => "\\\\")
            escaped_content = replace(escaped_content, "\"" => "\\\"")
            escaped_content = replace(escaped_content, "\n" => "\\n")
            escaped_content = replace(escaped_content, "\r" => "\\r")
            escaped_content = replace(escaped_content, "\t" => "\\t")
            call_msg = """{"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"scw_fixed_points","arguments":{"scw_content":"$escaped_content"}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 10
            @test haskey(resp, "result")
            content = resp["result"]["content"]
            @test length(content) >= 1
            text = content[1]["text"]

            # Same results as file_path test
            @test occursin("Consistent Scenarios (Fixed Points): 4", text)
            @test occursin("Total scenario space: 36", text)
            @test occursin("= 36 total", text)
        end

        # ── 6. tools/call — both file_path and scw_content (error) ──
        @testset "tools/call — both params error" begin
            scw_path = joinpath(SAMPLE_DIR, "CIB_global.scw")
            escaped_path = replace(scw_path, "\\" => "\\\\")
            call_msg = """{"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"scw_fixed_points","arguments":{"file_path":"$escaped_path","scw_content":"dummy"}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 11
            @test haskey(resp, "error")
        end

        # ── 7. tools/call — neither param (error) ──
        @testset "tools/call — no input error" begin
            call_msg = """{"jsonrpc":"2.0","id":12,"method":"tools/call","params":{"name":"scw_fixed_points","arguments":{}}}"""
            write(proc, mcp_frame(call_msg))
            flush(proc)

            resp = read_mcp_response(proc)
            @test resp["id"] == 12
            @test haskey(resp, "error")
        end

        # ── 8. tools/call with non-existent file ──
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

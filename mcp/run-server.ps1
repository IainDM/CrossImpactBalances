# PowerShell wrapper that launches the Node.js MCP server.
#
# The Node.js server handles the MCP protocol instantly and only
# spawns Julia when a tool is actually called, avoiding Julia's
# startup latency during the MCP handshake.
#
# Usage (Claude Code):
#   claude mcp add crossimpactbalances -- pwsh -NoProfile -File C:\path\to\CrossImpactBalances\mcp\run-server.ps1
#
# Usage (Claude Desktop / VS Code — claude_desktop_config.json):
#   {
#     "mcpServers": {
#       "crossimpactbalances": {
#         "command": "pwsh",
#         "args": ["-NoProfile", "-File", "C:\\path\\to\\CrossImpactBalances\\mcp\\run-server.ps1"]
#       }
#     }
#   }

$DIR = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
& node "$DIR\mcp\server.js"

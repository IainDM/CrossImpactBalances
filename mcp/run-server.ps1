# PowerShell wrapper that launches the MCP server.
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
& julia --startup-file=no --project="$DIR" "$DIR\mcp\server.jl"

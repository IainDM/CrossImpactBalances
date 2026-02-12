#!/bin/bash -l
# Wrapper that launches the Node.js MCP server under a login shell so
# that PATH includes both node and julia even when spawned by a
# non-interactive host (e.g. Claude Code / VS Code).
#
# The Node.js server handles the MCP protocol instantly and only
# spawns Julia when a tool is actually called, avoiding Julia's
# startup latency during the MCP handshake.
#
# Usage:
#   claude mcp add crossimpactbalances -- /absolute/path/to/CrossImpactBalances/mcp/run-server.sh

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec node "$DIR/mcp/server.js"

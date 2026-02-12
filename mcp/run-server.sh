#!/bin/bash -l
# Wrapper that launches the MCP server under a login shell so that
# PATH includes Julia even when spawned by a non-interactive host
# (e.g. Claude Code / VS Code).
#
# Usage:
#   claude mcp add crossimpactbalances -- /absolute/path/to/CrossImpactBalances/mcp/run-server.sh

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec julia --project="$DIR" "$DIR/mcp/server.jl"

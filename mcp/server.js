#!/usr/bin/env node
/**
 * Node.js MCP server wrapper for CrossImpactBalances.
 *
 * Handles the MCP protocol (initialize, tools/list) instantly in Node.js,
 * and spawns Julia only when a tool is actually called.  This avoids
 * Julia's startup latency blocking the MCP handshake.
 *
 * Zero dependencies — uses only Node.js built-ins.
 *
 * Usage:
 *   node mcp/server.js
 *
 * Claude Code:
 *   claude mcp add crossimpactbalances -- node /path/to/CrossImpactBalances/mcp/server.js
 *
 * Claude Desktop (claude_desktop_config.json):
 *   {
 *     "mcpServers": {
 *       "crossimpactbalances": {
 *         "command": "node",
 *         "args": ["/path/to/CrossImpactBalances/mcp/server.js"]
 *       }
 *     }
 *   }
 */
"use strict";

const { spawn } = require("child_process");
const path = require("path");

const PROJECT_DIR = path.resolve(__dirname, "..");
const ANALYSIS_SCRIPT = path.join(__dirname, "run_analysis.jl");

const SERVER_NAME = "crossimpactbalances-mcp";
const SERVER_VERSION = "0.1.0";
const PROTOCOL_VERSION = "2024-11-05";

// ── MCP transport (Content-Length framing over stdio) ──────────────────────

let inputBuffer = Buffer.alloc(0);
const pendingReads = [];

process.stdin.on("data", (chunk) => {
  inputBuffer = Buffer.concat([inputBuffer, chunk]);
  drainReads();
});

process.stdin.on("end", () => {
  for (const p of pendingReads) {
    p.reject(new Error("EOF"));
  }
  pendingReads.length = 0;
});

function tryParseMessage() {
  const sep = "\r\n\r\n";
  const headerEnd = inputBuffer.indexOf(sep);
  if (headerEnd === -1) return null;

  const headerStr = inputBuffer.slice(0, headerEnd).toString("utf-8");
  const match = headerStr.match(/Content-Length:\s*(\d+)/i);
  if (!match) {
    // Malformed header — skip it
    inputBuffer = inputBuffer.slice(headerEnd + sep.length);
    return null;
  }

  const contentLength = parseInt(match[1], 10);
  const bodyStart = headerEnd + sep.length;
  const bodyEnd = bodyStart + contentLength;

  if (inputBuffer.length < bodyEnd) return null; // incomplete body

  const body = inputBuffer.slice(bodyStart, bodyEnd).toString("utf-8");
  inputBuffer = inputBuffer.slice(bodyEnd);

  return JSON.parse(body);
}

function drainReads() {
  while (pendingReads.length > 0) {
    const msg = tryParseMessage();
    if (!msg) break;
    pendingReads.shift().resolve(msg);
  }
}

function readMessage() {
  const msg = tryParseMessage();
  if (msg) return Promise.resolve(msg);
  return new Promise((resolve, reject) => {
    pendingReads.push({ resolve, reject });
  });
}

function sendMessage(msg) {
  const body = JSON.stringify(msg);
  const header = `Content-Length: ${Buffer.byteLength(body, "utf-8")}\r\n\r\n`;
  process.stdout.write(header + body);
}

function log(text) {
  process.stderr.write(`${SERVER_NAME}: ${text}\n`);
}

// ── Protocol handlers ─────────────────────────────────────────────────────

function handleInitialize(id) {
  return {
    jsonrpc: "2.0",
    id,
    result: {
      protocolVersion: PROTOCOL_VERSION,
      capabilities: { tools: {} },
      serverInfo: { name: SERVER_NAME, version: SERVER_VERSION },
    },
  };
}

function handleToolsList(id) {
  return {
    jsonrpc: "2.0",
    id,
    result: {
      tools: [
        {
          name: "scw_fixed_points",
          description:
            "Analyze a ScenarioWizard (.scw) cross-impact balance file. " +
            "Finds all consistent scenarios (fixed points of the succession " +
            "operator) and their basins of attraction: the number of starting " +
            "combinations that converge to each fixed point.",
          inputSchema: {
            type: "object",
            properties: {
              file_path: {
                type: "string",
                description: "Absolute path to a ScenarioWizard .scw file",
              },
            },
            required: ["file_path"],
          },
        },
      ],
    },
  };
}

function handleToolsCall(id, params) {
  const name = (params && params.name) || "";
  const args = (params && params.arguments) || {};

  if (name === "scw_fixed_points") {
    return executeScwFixedPoints(id, args);
  }
  return Promise.resolve(makeError(id, -32601, `Unknown tool: ${name}`));
}

// ── Tool execution (spawns Julia) ─────────────────────────────────────────

function executeScwFixedPoints(id, args) {
  return new Promise((resolve) => {
    const filePath = args.file_path || "";
    if (!filePath) {
      resolve(makeError(id, -32602, "Missing required parameter: file_path"));
      return;
    }

    log(`spawning Julia for ${path.basename(filePath)}`);

    const proc = spawn(
      "julia",
      ["--startup-file=no", `--project=${PROJECT_DIR}`, ANALYSIS_SCRIPT, filePath],
      { stdio: ["ignore", "pipe", "pipe"] }
    );

    let stdout = "";
    let stderrBuf = "";

    proc.stdout.on("data", (d) => {
      stdout += d;
    });
    proc.stderr.on("data", (d) => {
      stderrBuf += d;
      process.stderr.write(d); // forward Julia's stderr to the MCP log
    });

    proc.on("error", (err) => {
      resolve(
        makeError(id, -32603, `Failed to start Julia: ${err.message}`)
      );
    });

    proc.on("close", (code) => {
      if (code === 0 && stdout) {
        log(`Julia finished successfully`);
        resolve({
          jsonrpc: "2.0",
          id,
          result: {
            content: [{ type: "text", text: stdout }],
          },
        });
      } else {
        const errMsg =
          stderrBuf.trim() || stdout.trim() || `Julia exited with code ${code}`;
        log(`Julia failed: ${errMsg}`);
        resolve(makeError(id, -32603, `Analysis failed: ${errMsg}`));
      }
    });
  });
}

function makeError(id, code, message) {
  return { jsonrpc: "2.0", id, error: { code, message } };
}

// ── Main loop ─────────────────────────────────────────────────────────────

async function main() {
  log(`starting (Node.js ${process.version})`);
  log("ready, waiting for messages");

  while (true) {
    let msg;
    try {
      msg = await readMessage();
    } catch (_) {
      break; // EOF or stdin closed
    }

    const method = msg.method || null;
    const id = msg.id !== undefined ? msg.id : null;
    const params = msg.params || {};

    log(`\u2190 ${method || "?"}${id !== null ? ` (id=${id})` : ""}`);

    let response = null;

    if (method === "initialize") {
      response = handleInitialize(id);
    } else if (method === "notifications/initialized") {
      // notification — no response
    } else if (method === "tools/list") {
      response = handleToolsList(id);
    } else if (method === "tools/call") {
      response = await handleToolsCall(id, params);
    } else if (id !== null) {
      response = makeError(id, -32601, `Method not found: ${method}`);
    }

    if (response) {
      sendMessage(response);
      log(`\u2192 response for ${method || "?"} (id=${id})`);
    }
  }

  log("shutting down (EOF)");
}

main().catch((err) => {
  log(`fatal error: ${err.message}`);
  process.exit(1);
});

#!/usr/bin/env node
/**
 * Node.js MCP server wrapper for CrossImpactBalances.
 *
 * Handles the MCP protocol (initialize, tools/list) instantly in Node.js.
 * On the first tool call, spawns a persistent Julia worker process that
 * loads CrossImpactBalances once and stays alive for subsequent calls.
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
const WORKER_SCRIPT = path.join(__dirname, "julia_worker.jl");

const SERVER_NAME = "crossimpactbalances-mcp";
const SERVER_VERSION = "0.1.0";
const PROTOCOL_VERSION = "2025-06-18";

// ── MCP transport (Content-Length framing over stdio) ──────────────────────

let inputBuffer = Buffer.alloc(0);
const pendingReads = [];

process.stdin.on("data", (chunk) => {
  log(`stdin: received ${chunk.length} bytes`);
  inputBuffer = Buffer.concat([inputBuffer, chunk]);
  drainReads();
});

process.stdin.on("error", (err) => {
  log(`stdin error: ${err.message}`);
});

process.stdin.on("end", () => {
  for (const p of pendingReads) {
    p.reject(new Error("EOF"));
  }
  pendingReads.length = 0;
});

function tryParseMessage() {
  // Support both \r\n\r\n (spec) and \n\n (some clients) as header delimiters
  const crlfIdx = inputBuffer.indexOf("\r\n\r\n");
  const lfIdx = inputBuffer.indexOf("\n\n");

  let headerEnd, sepLen;
  if (crlfIdx !== -1 && (lfIdx === -1 || crlfIdx <= lfIdx)) {
    headerEnd = crlfIdx;
    sepLen = 4; // \r\n\r\n
  } else if (lfIdx !== -1) {
    headerEnd = lfIdx;
    sepLen = 2; // \n\n
  } else {
    // No separator found — log buffer head for diagnostics
    if (inputBuffer.length > 0) {
      const preview = inputBuffer.slice(0, Math.min(80, inputBuffer.length));
      log(`buffer (no separator yet): ${JSON.stringify(preview.toString("utf-8"))}`);
    }
    return null;
  }

  const headerStr = inputBuffer.slice(0, headerEnd).toString("utf-8");
  const match = headerStr.match(/Content-Length:\s*(\d+)/i);
  if (!match) {
    // No Content-Length header — try parsing the whole buffer as bare JSON
    log(`no Content-Length header found in: ${JSON.stringify(headerStr)}`);
    inputBuffer = inputBuffer.slice(headerEnd + sepLen);
    return null;
  }

  const contentLength = parseInt(match[1], 10);
  const bodyStart = headerEnd + sepLen;
  const bodyEnd = bodyStart + contentLength;

  if (inputBuffer.length < bodyEnd) return null;

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

// ── Persistent Julia worker ───────────────────────────────────────────────

let workerProc = null;
let workerState = "idle"; // "idle" | "starting" | "ready"
let workerLineBuffer = "";

// Callbacks for the current in-flight operation
let onReady = null;   // { resolve, reject } — set while waiting for ready signal
let onResult = null;  // { resolve, reject } — set while waiting for analysis result

function processWorkerLines() {
  let idx;
  while ((idx = workerLineBuffer.indexOf("\n")) !== -1) {
    const line = workerLineBuffer.slice(0, idx).trim();
    workerLineBuffer = workerLineBuffer.slice(idx + 1);
    if (!line) continue;

    let parsed;
    try {
      parsed = JSON.parse(line);
    } catch (_) {
      continue; // ignore non-JSON lines
    }

    if (workerState === "starting" && parsed.status === "ready") {
      workerState = "ready";
      log("Julia worker ready");
      if (onReady) {
        onReady.resolve();
        onReady = null;
      }
    } else if (workerState === "ready" && onResult) {
      const cb = onResult;
      onResult = null;
      cb.resolve(parsed);
    }
  }
}

function cleanupWorker() {
  workerProc = null;
  workerState = "idle";
  workerLineBuffer = "";

  const err = new Error("Julia worker exited unexpectedly");
  if (onReady) {
    onReady.reject(err);
    onReady = null;
  }
  if (onResult) {
    onResult.reject(err);
    onResult = null;
  }
}

function ensureWorker() {
  if (workerState === "ready") return Promise.resolve();

  if (workerState === "starting") {
    return new Promise((resolve, reject) => {
      onReady = { resolve, reject };
    });
  }

  // Spawn a new worker
  return new Promise((resolve, reject) => {
    workerState = "starting";
    onReady = { resolve, reject };

    log("spawning persistent Julia worker");

    workerProc = spawn(
      "julia",
      ["--startup-file=no", `--project=${PROJECT_DIR}`, WORKER_SCRIPT],
      { stdio: ["pipe", "pipe", "pipe"] }
    );

    workerProc.stderr.on("data", (d) => {
      process.stderr.write(d); // forward Julia's stderr to MCP log
    });

    workerProc.stdout.on("data", (chunk) => {
      workerLineBuffer += chunk.toString();
      processWorkerLines();
    });

    workerProc.on("error", (err) => {
      log(`Julia worker spawn error: ${err.message}`);
      cleanupWorker();
    });

    workerProc.on("close", (code) => {
      log(`Julia worker exited (code ${code})`);
      cleanupWorker();
    });
  });
}

function sendToWorker(filePath) {
  return new Promise((resolve, reject) => {
    onResult = { resolve, reject };
    workerProc.stdin.write(filePath + "\n");
  });
}

// Kill worker on exit
process.on("exit", () => {
  if (workerProc) {
    workerProc.kill();
  }
});

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

async function handleToolsCall(id, params) {
  const name = (params && params.name) || "";
  const args = (params && params.arguments) || {};

  if (name === "scw_fixed_points") {
    return executeScwFixedPoints(id, args);
  }
  return makeError(id, -32601, `Unknown tool: ${name}`);
}

// ── Tool execution (uses persistent Julia worker) ─────────────────────────

async function executeScwFixedPoints(id, args) {
  const filePath = args.file_path || "";
  if (!filePath) {
    return makeError(id, -32602, "Missing required parameter: file_path");
  }

  try {
    await ensureWorker();

    log(`analyzing ${path.basename(filePath)}`);
    const result = await sendToWorker(filePath);

    if (result.ok) {
      log("analysis complete");
      return {
        jsonrpc: "2.0",
        id,
        result: {
          content: [{ type: "text", text: result.text }],
        },
      };
    } else {
      return makeError(id, -32603, `Analysis failed: ${result.error}`);
    }
  } catch (err) {
    return makeError(id, -32603, `Julia worker error: ${err.message}`);
  }
}

function makeError(id, code, message) {
  return { jsonrpc: "2.0", id, error: { code, message } };
}

// ── Main loop ─────────────────────────────────────────────────────────────

async function main() {
  log(`starting (Node.js ${process.version})`);
  log(`stdin: isTTY=${process.stdin.isTTY}, readable=${process.stdin.readable}`);
  process.stdin.resume(); // ensure flowing mode on Windows
  log("ready, waiting for messages");

  while (true) {
    let msg;
    try {
      msg = await readMessage();
    } catch (_) {
      break;
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
      const body = JSON.stringify(response);
      log(`stdout: sending ${Buffer.byteLength(body)} bytes for ${method || "?"} (id=${id})`);
      sendMessage(response);
      log(`\u2192 response sent for ${method || "?"} (id=${id})`);
    }
  }

  log("shutting down (EOF)");
}

main().catch((err) => {
  log(`fatal error: ${err.message}`);
  process.exit(1);
});

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

// ── MCP transport (auto-detects Content-Length framing vs bare NDJSON) ─────

let inputBuffer = Buffer.alloc(0);
const pendingReads = [];
let useContentLength = null; // auto-detected from first incoming message

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
  if (inputBuffer.length === 0) return null;

  // Skip leading whitespace / empty lines
  let start = 0;
  while (start < inputBuffer.length) {
    const b = inputBuffer[start];
    if (b !== 0x20 && b !== 0x09 && b !== 0x0a && b !== 0x0d) break;
    start++;
  }
  if (start > 0) inputBuffer = inputBuffer.slice(start);
  if (inputBuffer.length === 0) return null;

  // Auto-detect transport: '{' means bare JSON, otherwise Content-Length
  if (inputBuffer[0] === 0x7b) {
    // ── Bare JSON (newline-delimited) ──
    const nlIdx = inputBuffer.indexOf(0x0a);
    if (nlIdx === -1) return null; // wait for complete line

    const line = inputBuffer.slice(0, nlIdx).toString("utf-8").trim();
    inputBuffer = inputBuffer.slice(nlIdx + 1);

    if (!line) return null;

    if (useContentLength === null) {
      useContentLength = false;
      log("transport: bare JSON (NDJSON) detected");
    }
    return JSON.parse(line);
  }

  // ── Content-Length framed ──
  const crlfIdx = inputBuffer.indexOf("\r\n\r\n");
  const lfIdx = inputBuffer.indexOf("\n\n");

  let headerEnd, sepLen;
  if (crlfIdx !== -1 && (lfIdx === -1 || crlfIdx <= lfIdx)) {
    headerEnd = crlfIdx;
    sepLen = 4;
  } else if (lfIdx !== -1) {
    headerEnd = lfIdx;
    sepLen = 2;
  } else {
    return null; // wait for complete header
  }

  const headerStr = inputBuffer.slice(0, headerEnd).toString("utf-8");
  const match = headerStr.match(/Content-Length:\s*(\d+)/i);
  if (!match) {
    inputBuffer = inputBuffer.slice(headerEnd + sepLen);
    return null;
  }

  const contentLength = parseInt(match[1], 10);
  const bodyStart = headerEnd + sepLen;
  const bodyEnd = bodyStart + contentLength;

  if (inputBuffer.length < bodyEnd) return null;

  const body = inputBuffer.slice(bodyStart, bodyEnd).toString("utf-8");
  inputBuffer = inputBuffer.slice(bodyEnd);

  if (useContentLength === null) {
    useContentLength = true;
    log("transport: Content-Length framing detected");
  }
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
  if (useContentLength) {
    const header = `Content-Length: ${Buffer.byteLength(body, "utf-8")}\r\n\r\n`;
    process.stdout.write(header + body);
  } else {
    process.stdout.write(body + "\n");
  }
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

function sendToWorker(op, filePath, scenario) {
  return new Promise((resolve, reject) => {
    onResult = { resolve, reject };
    let cmd = `${op}\t${filePath}`;
    if (scenario) {
      for (const [desc, variant] of Object.entries(scenario)) {
        cmd += `\t${desc}=${variant}`;
      }
    }
    workerProc.stdin.write(cmd + "\n");
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
            "Find all consistent scenarios (fixed points) in a " +
            "ScenarioWizard (.scw) cross-impact balance file. Lists the " +
            "scenarios that are self-consistent under the succession operator. " +
            "Fast — use scw_basins instead if you also need basin sizes.",
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
        {
          name: "scw_basins",
          description:
            "Exhaustive basin-of-attraction analysis of a ScenarioWizard " +
            "(.scw) file. Finds all consistent scenarios AND counts how " +
            "many starting combinations converge to each one. Slower than " +
            "scw_fixed_points but gives the full convergence picture.",
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
        {
          name: "scw_succession",
          description:
            "Trace the succession trajectory from a given starting scenario. " +
            "Shows step-by-step how the scenario evolves under the CIB " +
            "succession operator until it converges to a fixed point or " +
            "enters a cycle. Use this to explore 'what if we start from X?'.",
          inputSchema: {
            type: "object",
            properties: {
              file_path: {
                type: "string",
                description: "Absolute path to a ScenarioWizard .scw file",
              },
              scenario: {
                type: "object",
                description:
                  "The starting scenario as a mapping of descriptor names " +
                  "to variant names. Every descriptor must be included. " +
                  'Example: {"Economy": "Growth", "Climate": "Warming"}',
                additionalProperties: { type: "string" },
              },
            },
            required: ["file_path", "scenario"],
          },
        },
        {
          name: "scw_impact_balance",
          description:
            "Compute the impact balance for a given scenario, showing the " +
            "support score for every variant of every descriptor. Reveals " +
            "why a scenario is or isn't self-consistent: which variants " +
            "are favoured and which would change under succession.",
          inputSchema: {
            type: "object",
            properties: {
              file_path: {
                type: "string",
                description: "Absolute path to a ScenarioWizard .scw file",
              },
              scenario: {
                type: "object",
                description:
                  "The scenario to analyze as a mapping of descriptor names " +
                  "to variant names. Every descriptor must be included. " +
                  'Example: {"Economy": "Growth", "Climate": "Warming"}',
                additionalProperties: { type: "string" },
              },
            },
            required: ["file_path", "scenario"],
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
    return executeWorkerTool(id, "fixed_points", args);
  }
  if (name === "scw_basins") {
    return executeWorkerTool(id, "basins", args);
  }
  if (name === "scw_succession") {
    return executeWorkerTool(id, "succession", args);
  }
  if (name === "scw_impact_balance") {
    return executeWorkerTool(id, "impact_balance", args);
  }
  return makeError(id, -32601, `Unknown tool: ${name}`);
}

// ── Tool execution (uses persistent Julia worker) ─────────────────────────

async function executeWorkerTool(id, op, args) {
  const filePath = args.file_path || "";
  if (!filePath) {
    return makeError(id, -32602, "Missing required parameter: file_path");
  }

  const scenario = args.scenario || null;
  if ((op === "succession" || op === "impact_balance") && !scenario) {
    return makeError(id, -32602, "Missing required parameter: scenario");
  }

  try {
    await ensureWorker();

    log(`${op} on ${path.basename(filePath)}`);
    const result = await sendToWorker(op, filePath, scenario);

    if (result.ok) {
      log(`${op} complete`);
      return {
        jsonrpc: "2.0",
        id,
        result: {
          content: [{ type: "text", text: result.text }],
        },
      };
    } else {
      return makeError(id, -32603, `${result.error}`);
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

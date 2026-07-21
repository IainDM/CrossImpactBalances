"""
    CIBApp

A small local desktop application for Cross-Impact Balance scenario analysis.

`julia_main()` starts a localhost web server and opens the default browser to
a single-page UI: pick a ScenarioWizard `.scw` file, then run *Find Consistent
Scenarios* or *Find Basins*. Results are shown as a table (with basin sizes for
the basin analysis) and the full basin analysis can be exported as CSV.

The server has no external dependencies — it is a minimal HTTP/1.1 responder
built on the standard-library `Sockets`, talking only to the bundled page — so
the whole application compiles into a single self-contained executable.
"""
module CIBApp

using Sockets
using CrossImpactBalances

# ─── Server state (single local user) ───────────────────────────────────────
mutable struct AppState
    last_csv::String     # CSV of the most recent basin analysis, for export
    shutdown::Bool
end
AppState() = AppState("", false)

# ─── Minimal JSON encoding (we control every value that is encoded) ──────────
function json_escape(s::AbstractString)
    io = IOBuffer()
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif c < ' '
            print(io, "\\u", lpad(string(UInt16(c), base = 16), 4, '0'))
        else
            print(io, c)
        end
    end
    return String(take!(io))
end

to_json(s::AbstractString) = string('"', json_escape(s), '"')
to_json(b::Bool) = b ? "true" : "false"
to_json(n::Integer) = string(n)
to_json(x::AbstractFloat) = isfinite(x) ? string(x) : "null"
to_json(::Nothing) = "null"
to_json(v::AbstractVector) = string('[', join(map(to_json, v), ','), ']')
function to_json(d::AbstractDict)
    parts = String[string(to_json(string(k)), ':', to_json(v)) for (k, v) in d]
    return string('{', join(parts, ','), '}')
end

# ─── Minimal CSV field quoting ───────────────────────────────────────────────
function csv_field(s::AbstractString)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return string('"', replace(s, '"' => "\"\""), '"')
    end
    return String(s)
end
csv_field(x) = csv_field(string(x))
csv_row(fields) = join(map(csv_field, fields), ',')

# ─── Analysis ────────────────────────────────────────────────────────────────
# Variant names of a scenario (0-based indices) as human-readable strings.
function variant_names(cib, u)
    return String[cib.variants[cib.descriptors[i]][u[i] + 1] for i in 1:cib.ndesc]
end

# Load a CIB from raw .scw text without triggering any automatic search
# (an empty, non-nothing kernel suppresses load_scw's own find_consistent).
function load_from_text(scwtext::AbstractString)
    path = tempname() * ".scw"
    write(path, scwtext)
    try
        return load_scw(path; kernel = Vector{Vector{Int}}())
    finally
        rm(path; force = true)
    end
end

function analyze_consistent(cib)
    kern = find_consistent(cib; exhaustive = true)
    total = max_signature(cib) + 1
    scenarios = [Dict("signature" => signature(cib, u),
                      "variants" => variant_names(cib, u)) for u in kern]
    return Dict("mode" => "consistent",
                "descriptors" => cib.descriptors,
                "total" => total,
                "count" => length(kern),
                "scenarios" => scenarios)
end

function analyze_basins(cib, state::AppState)
    fps, sizes, cycles = find_basins(cib)
    total = max_signature(cib) + 1
    order = sortperm(sizes; rev = true)      # largest basin first

    scenarios = Dict{String,Any}[]
    for k in order
        u = fps[k]
        push!(scenarios, Dict("signature" => signature(cib, u),
                              "variants" => variant_names(cib, u),
                              "basin_size" => sizes[k],
                              "basin_pct" => round(100 * sizes[k] / total; digits = 4)))
    end

    # Build the export CSV and stash it for /export.csv.
    io = IOBuffer()
    println(io, "# Cross-Impact Balances basin analysis")
    println(io, "# total scenarios,", total)
    println(io, "# consistent scenarios,", length(fps))
    println(io, "# into-cycle starts,", cycles)
    println(io, csv_row(vcat(["rank", "signature"], cib.descriptors,
                             ["basin_size", "basin_fraction_pct"])))
    for (rank, k) in enumerate(order)
        u = fps[k]
        pct = round(100 * sizes[k] / total; digits = 4)
        # Any[] avoids Int→Float promotion so basin_size stays an integer.
        println(io, csv_row(vcat(Any[rank, signature(cib, u)], variant_names(cib, u),
                                 Any[sizes[k], pct])))
    end
    state.last_csv = String(take!(io))

    return Dict("mode" => "basins",
                "descriptors" => cib.descriptors,
                "total" => total,
                "count" => length(fps),
                "cycles" => cycles,
                "covered" => sum(sizes; init = 0),
                "scenarios" => scenarios)
end

# ─── Minimal HTTP/1.1 ────────────────────────────────────────────────────────
# Read one request (request line, headers, and a Content-Length body). Returns
# (method, target, body) or nothing on a closed/blank connection.
function read_request(sock)
    reqline = try
        readline(sock)
    catch
        return nothing
    end
    isempty(reqline) && return nothing
    parts = split(reqline)
    length(parts) < 2 && return nothing
    method = String(parts[1])
    target = String(parts[2])

    clen = 0
    while true
        line = readline(sock)
        isempty(line) && break                 # blank line ends the headers
        idx = findfirst(':', line)
        idx === nothing && continue
        if lowercase(strip(line[1:idx-1])) == "content-length"
            clen = something(tryparse(Int, strip(line[idx+1:end])), 0)
        end
    end

    body = ""
    if clen > 0
        body = String(read(sock, clen))
    end
    return (method, target, body)
end

const STATUS_TEXT = Dict(200 => "OK", 204 => "No Content", 400 => "Bad Request",
                         404 => "Not Found", 500 => "Internal Server Error")

function respond(sock, status::Int, ctype::String, body;
                 extra::Vector{Pair{String,String}} = Pair{String,String}[])
    data = body isa Vector{UInt8} ? body : Vector{UInt8}(codeunits(String(body)))
    io = IOBuffer()
    print(io, "HTTP/1.1 ", status, ' ', get(STATUS_TEXT, status, "OK"), "\r\n")
    print(io, "Content-Type: ", ctype, "\r\n")
    print(io, "Content-Length: ", length(data), "\r\n")
    for (k, v) in extra
        print(io, k, ": ", v, "\r\n")
    end
    print(io, "Connection: close\r\n\r\n")
    write(sock, take!(io))
    isempty(data) || write(sock, data)
    return nothing
end

# Split "/analyze?mode=basins" into ("/analyze", Dict("mode"=>"basins")).
function split_target(target)
    q = findfirst('?', target)
    q === nothing && return (target, Dict{String,String}())
    path = target[1:q-1]
    params = Dict{String,String}()
    for kv in split(target[q+1:end], '&')
        isempty(kv) && continue
        eq = findfirst('=', kv)
        if eq === nothing
            params[String(kv)] = ""
        else
            params[String(kv[1:eq-1])] = String(kv[eq+1:end])
        end
    end
    return (path, params)
end

function handle(sock, state::AppState)
    req = read_request(sock)
    req === nothing && return
    method, target, body = req
    path, params = split_target(target)

    if method == "GET" && (path == "/" || path == "/index.html")
        respond(sock, 200, "text/html; charset=utf-8", PAGE)
    elseif method == "GET" && path == "/favicon.ico"
        respond(sock, 204, "text/plain", "")
    elseif method == "GET" && path == "/health"
        respond(sock, 200, "text/plain", "ok")
    elseif method == "GET" && path == "/quit"
        state.shutdown = true
        respond(sock, 200, "text/plain", "Shutting down. You can close this tab.")
    elseif method == "GET" && path == "/export.csv"
        if isempty(state.last_csv)
            respond(sock, 400, "text/plain", "Run a basin analysis first.")
        else
            respond(sock, 200, "text/csv; charset=utf-8", state.last_csv;
                    extra = ["Content-Disposition" =>
                             "attachment; filename=\"basin_analysis.csv\""])
        end
    elseif method == "POST" && path == "/analyze"
        mode = get(params, "mode", "consistent")
        try
            cib = load_from_text(body)
            result = mode == "basins" ? analyze_basins(cib, state) :
                     analyze_consistent(cib)
            respond(sock, 200, "application/json; charset=utf-8", to_json(result))
        catch e
            msg = sprint(showerror, e)
            respond(sock, 400, "application/json; charset=utf-8",
                    to_json(Dict("error" => msg)))
        end
    else
        respond(sock, 404, "text/plain", "Not found")
    end
    return nothing
end

# ─── Server loop + entry point ───────────────────────────────────────────────
function open_browser(url)
    if Sys.iswindows()
        run(`cmd /c start "" $url`)
    elseif Sys.isapple()
        run(`open $url`)
    else
        run(`xdg-open $url`)
    end
    return nothing
end

function serve(state::AppState, server)
    while !state.shutdown
        sock = try
            accept(server)
        catch
            break                               # listener closed
        end
        try
            handle(sock, state)
        catch e
            try
                respond(sock, 500, "text/plain", "Internal error: " * sprint(showerror, e))
            catch
            end
        finally
            close(sock)
        end
    end
    try
        close(server)
    catch
    end
    return nothing
end

"""
    julia_main() -> Cint

Application entry point (also the PackageCompiler target). Binds a localhost
port, opens the browser, and serves the UI until the browser's Quit button is
used or the process is interrupted.
"""
function julia_main()::Cint
    server = nothing
    port = 0
    for p in 8071:8099
        try
            server = listen(Sockets.localhost, p)
            port = p
            break
        catch
            continue
        end
    end
    if server === nothing
        println(stderr, "CIBApp: could not bind a local port in 8071–8099.")
        return 1
    end

    url = "http://127.0.0.1:$port/"
    println("""
    ============================================================
      Cross-Impact Balances — desktop app
      Open in your browser:  $url
      Threads available:     $(Threads.nthreads())
      Leave this window open while you work; use the Quit
      button in the page (or close this window) to stop.
    ============================================================
    """)
    flush(stdout)

    try
        open_browser(url)
    catch
        println("(Could not open the browser automatically — open $url manually.)")
    end

    serve(AppState(), server)
    println("CIBApp stopped.")
    return 0
end

# ─── Embedded single-page UI ─────────────────────────────────────────────────
const PAGE = raw"""
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cross-Impact Balances</title>
<style>
  :root { --bg:#f6f7f9; --card:#fff; --line:#e3e6ea; --ink:#1c2530;
          --muted:#6b7683; --accent:#2a6df4; --accent2:#0f9d58; }
  * { box-sizing: border-box; }
  body { margin:0; background:var(--bg); color:var(--ink);
         font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; }
  header { background:var(--card); border-bottom:1px solid var(--line);
           padding:18px 24px; }
  header h1 { margin:0; font-size:19px; }
  header p { margin:4px 0 0; color:var(--muted); font-size:13px; }
  main { max-width:1100px; margin:22px auto; padding:0 24px; }
  .card { background:var(--card); border:1px solid var(--line); border-radius:10px;
          padding:18px 20px; margin-bottom:18px; }
  .row { display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
  .filebtn { position:relative; overflow:hidden; display:inline-block; }
  input[type=file] { position:absolute; left:0; top:0; opacity:0; width:100%;
                     height:100%; cursor:pointer; }
  button { font:inherit; border:1px solid var(--line); background:#fff;
           padding:9px 16px; border-radius:8px; cursor:pointer; color:var(--ink); }
  button:hover:not(:disabled) { border-color:#c3c9d1; }
  button:disabled { opacity:.5; cursor:not-allowed; }
  button.primary { background:var(--accent); border-color:var(--accent); color:#fff; }
  button.green   { background:var(--accent2); border-color:var(--accent2); color:#fff; }
  .fname { color:var(--muted); font-size:13px; }
  .muted { color:var(--muted); }
  #status { margin-top:12px; min-height:20px; font-size:14px; }
  .spinner { display:inline-block; width:14px; height:14px; border:2px solid #c9d2df;
             border-top-color:var(--accent); border-radius:50%;
             animation:spin .8s linear infinite; vertical-align:-2px; margin-right:6px; }
  @keyframes spin { to { transform:rotate(360deg); } }
  .summary { margin:2px 0 14px; font-size:14px; }
  .summary b { font-size:16px; }
  .tablewrap { overflow-x:auto; border:1px solid var(--line); border-radius:8px; }
  table { border-collapse:collapse; width:100%; font-size:13px; white-space:nowrap; }
  th, td { padding:7px 11px; border-bottom:1px solid var(--line); text-align:left; }
  th { background:#f0f2f5; position:sticky; top:0; font-weight:600; }
  td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; }
  tr:last-child td { border-bottom:none; }
  .bar { height:8px; background:var(--accent); border-radius:4px; }
  .err { color:#b3261e; }
  .toolbar { display:flex; gap:10px; align-items:center; margin-bottom:12px;
             flex-wrap:wrap; }
  footer { text-align:center; color:var(--muted); font-size:12px; padding:10px 0 30px; }
</style>
</head>
<body>
<header>
  <h1>Cross-Impact Balances</h1>
  <p>Load a ScenarioWizard <code>.scw</code> file, then find its consistent
     scenarios or run a full basin-of-attraction analysis.</p>
</header>
<main>
  <div class="card">
    <div class="row">
      <span class="filebtn"><button type="button">Browse…</button>
        <input type="file" id="file" accept=".scw,.txt"></span>
      <span class="fname" id="fname">No file selected</span>
    </div>
    <div class="row" style="margin-top:14px;">
      <button class="primary" id="btnConsistent" disabled>Find Consistent Scenarios</button>
      <button class="green" id="btnBasins" disabled>Find Basins</button>
      <span style="flex:1"></span>
      <button id="btnQuit" title="Stop the app">Quit</button>
    </div>
    <div id="status"></div>
  </div>
  <div class="card" id="results" style="display:none;">
    <div class="toolbar">
      <div class="summary" id="summary" style="flex:1"></div>
      <a id="exportLink" href="/export.csv" style="display:none;">
        <button class="green">Export basin CSV</button></a>
    </div>
    <div class="tablewrap"><table id="table"></table></div>
  </div>
  <footer>Cross-Impact Balances desktop app · runs locally on your machine</footer>
</main>
<script>
let scwText = null;
const $ = id => document.getElementById(id);

$("file").addEventListener("change", e => {
  const f = e.target.files[0];
  if (!f) return;
  const reader = new FileReader();
  reader.onload = () => {
    scwText = reader.result;
    $("fname").textContent = f.name + " (" + f.size.toLocaleString() + " bytes)";
    $("btnConsistent").disabled = false;
    $("btnBasins").disabled = false;
    $("status").textContent = "";
  };
  reader.readAsText(f);
});

function setBusy(msg) {
  $("status").innerHTML = '<span class="spinner"></span>' + msg;
  $("btnConsistent").disabled = true;
  $("btnBasins").disabled = true;
}
function clearBusy() {
  $("status").textContent = "";
  $("btnConsistent").disabled = false;
  $("btnBasins").disabled = false;
}

async function analyze(mode) {
  if (!scwText) return;
  setBusy(mode === "basins"
      ? "Running full basin analysis… (this can take a moment on large files)"
      : "Searching for consistent scenarios…");
  try {
    const t0 = performance.now();
    const res = await fetch("/analyze?mode=" + mode, { method:"POST", body: scwText });
    const data = await res.json();
    const secs = ((performance.now() - t0) / 1000).toFixed(3);
    if (data.error) { showError(data.error); return; }
    render(data, secs);
  } catch (err) {
    showError(err.message || String(err));
  } finally {
    clearBusy();
  }
}

function showError(msg) {
  $("results").style.display = "none";
  $("status").innerHTML = '<span class="err">Error: ' + escapeHtml(msg) + '</span>';
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

function render(data, secs) {
  const isBasins = data.mode === "basins";
  const maxBasin = isBasins && data.scenarios.length
      ? Math.max(...data.scenarios.map(s => s.basin_size)) : 0;

  let sum;
  if (isBasins) {
    const cov = (100 * data.covered / data.total).toFixed(2);
    sum = `<b>${data.count}</b> consistent scenario(s) across a space of `
        + `${data.total.toLocaleString()}. They attract <b>${cov}%</b> of all `
        + `scenarios; ${data.cycles.toLocaleString()} start(s) fall into `
        + `non-fixed-point cycles. <span class="muted">(${secs}s)</span>`;
    $("exportLink").style.display = "inline-block";
  } else {
    sum = `<b>${data.count}</b> consistent scenario(s) out of `
        + `${data.total.toLocaleString()}. `
        + `<span class="muted">(${secs}s)</span>`;
    $("exportLink").style.display = "none";
  }
  $("summary").innerHTML = sum;

  const cols = ["#", "Signature", ...data.descriptors];
  if (isBasins) cols.push("Basin size", "Basin %", "");
  let html = "<thead><tr>";
  cols.forEach((c, i) => {
    const numeric = c === "Signature" || c === "Basin size" || c === "Basin %";
    html += `<th class="${numeric ? 'num' : ''}">${escapeHtml(c)}</th>`;
  });
  html += "</tr></thead><tbody>";

  data.scenarios.forEach((s, i) => {
    html += "<tr>";
    html += `<td class="num">${i + 1}</td>`;
    html += `<td class="num">${s.signature.toLocaleString()}</td>`;
    s.variants.forEach(v => html += `<td>${escapeHtml(v)}</td>`);
    if (isBasins) {
      html += `<td class="num">${s.basin_size.toLocaleString()}</td>`;
      html += `<td class="num">${s.basin_pct}</td>`;
      const w = maxBasin ? Math.max(2, Math.round(70 * s.basin_size / maxBasin)) : 0;
      html += `<td><div class="bar" style="width:${w}px"></div></td>`;
    }
    html += "</tr>";
  });
  html += "</tbody>";
  $("table").innerHTML = html;
  $("results").style.display = "block";
}

$("btnConsistent").addEventListener("click", () => analyze("consistent"));
$("btnBasins").addEventListener("click", () => analyze("basins"));
$("btnQuit").addEventListener("click", async () => {
  try { await fetch("/quit"); } catch (e) {}
  document.body.innerHTML =
    '<main><div class="card">The app has stopped. You can close this tab.</div></main>';
});
</script>
</body>
</html>
"""

end # module

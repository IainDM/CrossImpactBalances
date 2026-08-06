"""
    CIBApp

A small local desktop application for Cross-Impact Balance scenario analysis.

`julia_main()` starts a localhost web server and opens the default browser to
a single-page UI: pick a ScenarioWizard `.scw` file, then run *Find Consistent
Scenarios*, *Exact Basins*, *Estimate Basin Shares* (sampling with confidence
intervals, for spaces too large for any exact method — the consistent
scenarios themselves are still found exactly first), or *Structure* (the
influence map: independent islands, dial descriptors, and — when the model
decomposes into small enough islands — exact composed basins on the spot).
Results are shown as a table; basin analyses and estimates export as CSV. A
model too large for the exact table analysis gets guidance and a one-click
path to the estimator instead of an error.

The server has no external dependencies — it is a minimal HTTP/1.1 responder
built on the standard-library `Sockets`, talking only to the bundled page — so
the whole application compiles into a single self-contained executable.
"""
module CIBApp

using Sockets
using CrossImpactBalances

# ─── Server state (single local user) ───────────────────────────────────────
mutable struct AppState
    last_csv::String       # CSV of the most recent basin/estimate analysis, for export
    last_csv_name::String  # download filename for it (exact vs estimated differ)
    shutdown::Bool
end
AppState() = AppState("", "basin_analysis.csv", false)

# Counts stay JSON numbers while a double can hold them exactly; past 2^53
# they cross as strings so the browser cannot mangle their last digits.
# (Int128-safe: structure/estimate results carry counts far beyond Int64.)
json_count(n::Integer) = Int128(n) <= Int128(2)^53 ? Int(n) : string(n)

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
    return String[cib.variants[cib.descriptors[i]][u[i] + 1] for i in 1:cib.numberOfDescriptors]
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
    compute = @elapsed (kern = find_consistent(cib))
    total = max_signature(cib) + 1
    scenarios = [Dict("signature" => signature(cib, u),
                      "variants" => variant_names(cib, u)) for u in kern]
    return Dict("mode" => "consistent",
                "descriptors" => cib.descriptors,
                # Past 2^53 a JSON number silently loses its last digits in the
                # browser (JS numbers are doubles); a string crosses intact and
                # renders unchanged. Real matrices with >9e15 scenarios exist.
                "total" => total > 2^53 ? string(total) : total,
                "count" => length(kern),
                "compute_s" => round(compute; digits = 3),
                "threads" => Threads.nthreads(),
                "scenarios" => scenarios)
end

function analyze_basins(cib, state::AppState)
    compute = @elapsed ((fps, sizes, cycles) = find_basins(cib))
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
    state.last_csv_name = "basin_analysis.csv"

    return Dict("mode" => "basins",
                "descriptors" => cib.descriptors,
                "total" => total,
                "count" => length(fps),
                "cycles" => cycles,
                "covered" => sum(sizes; init = 0),
                "compute_s" => round(compute; digits = 3),
                "threads" => Threads.nthreads(),
                "scenarios" => scenarios)
end

# Basin SHARES by sampling, for models whose space defeats every exact method.
# The kernel is still found exactly first (estimate_basins runs find_consistent
# itself when the loaded model carries no kernel), so every consistent scenario
# appears in the report — a never-hit one shows its upper bound, not a false 0.
function analyze_estimate(cib, state::AppState; samples::Int = 1_000_000)
    compute = @elapsed (est = estimate_basins(cib; samples = samples))
    order = sortperm(est.hits; rev = true)     # biggest estimated share first

    scenarios = Dict{String,Any}[]
    for k in order
        push!(scenarios, Dict(
            "variants" => variant_names(cib, est.fixedPoints[k]),
            "hits" => est.hits[k],
            "share_pct" => round(100 * est.shares[k]; sigdigits = 4),
            "ci_lo_pct" => round(100 * est.ciLow[k]; sigdigits = 4),
            "ci_hi_pct" => round(100 * est.ciHigh[k]; sigdigits = 4),
            "est_size" => round(est.sizeEstimates[k]; sigdigits = 3)))
    end

    io = IOBuffer()
    println(io, "# Cross-Impact Balances basin-share estimate (exact kernel, sampled shares)")
    println(io, "# total scenarios,", est.scenarioCount)
    println(io, "# samples,", est.samples)
    println(io, "# seed,0x", string(est.seed, base = 16))
    println(io, "# confidence,", est.confidence)
    println(io, "# cycle share pct,", round(100 * est.cycleShare; sigdigits = 4),
            ",ci_low,", round(100 * est.cycleCiLow; sigdigits = 4),
            ",ci_high,", round(100 * est.cycleCiHigh; sigdigits = 4))
    println(io, csv_row(vcat(["rank"], cib.descriptors,
                             ["hits", "share_pct", "ci_low_pct", "ci_high_pct",
                              "estimated_size"])))
    for (rank, k) in enumerate(order)
        println(io, csv_row(vcat(Any[rank], variant_names(cib, est.fixedPoints[k]),
                                 Any[est.hits[k],
                                     round(100 * est.shares[k]; sigdigits = 4),
                                     round(100 * est.ciLow[k]; sigdigits = 4),
                                     round(100 * est.ciHigh[k]; sigdigits = 4),
                                     round(est.sizeEstimates[k]; sigdigits = 3)])))
    end
    state.last_csv = String(take!(io))
    state.last_csv_name = "basin_share_estimate.csv"

    return Dict("mode" => "estimate",
                "descriptors" => cib.descriptors,
                "total" => json_count(est.scenarioCount),
                "samples" => est.samples,
                "seed" => "0x" * string(est.seed, base = 16),
                "confidence" => est.confidence,
                "count" => length(scenarios),
                "cycle_hits" => est.cycleHits,
                "cycle_share_pct" => round(100 * est.cycleShare; sigdigits = 4),
                "cycle_ci_lo_pct" => round(100 * est.cycleCiLow; sigdigits = 4),
                "cycle_ci_hi_pct" => round(100 * est.cycleCiHigh; sigdigits = 4),
                "compute_s" => round(compute; digits = 3),
                "threads" => Threads.nthreads(),
                "scenarios" => scenarios)
end

# The influence map: who really has a say over whom, the independent islands,
# and the dial / jump-and-stay descriptors. When the model decomposes into
# islands that are individually small enough, the exact composed basins are
# computed on the spot (product_basins) — exact answers for spaces far past
# what the flat-table analysis could touch.
function analyze_structure(cib)
    compute = @elapsed (st = influence_structure(cib))
    islandScenarios(component) = prod(Int128.(cib.numberOfVariants[component]))

    islands = [Dict("descriptors" => st.descriptors[component],
                    "scenarios" => json_count(islandScenarios(component)))
               for component in st.components]
    dials = [Dict("descriptor" => st.descriptors[d],
                  "variants" => cib.variants[st.descriptors[d]]) for d in st.tieFrozen]
    forced = [Dict("descriptor" => st.descriptors[d],
                   "variants" => [cib.variants[st.descriptors[d]][v + 1] for v in maximizers])
              for (d, maximizers) in st.forced]

    result = Dict{String,Any}(
        "mode" => "structure",
        "descriptors" => cib.descriptors,
        "total" => json_count(scenario_count(cib)),
        "influences" => count(st.activeEdges),
        "islands" => islands,
        "dials" => dials,
        "forced" => forced,
        "threads" => Threads.nthreads())

    if length(st.components) > 1
        largestIsland = maximum(islandScenarios(c) for c in st.components)
        if largestIsland <= 50_000_000
            # Every island is a sub-second table job: compose the exact answer.
            # (product_basins refuses degenerate models whose island kernels
            # combine to astronomically many scenarios — that refusal is
            # guidance here, not a failure.)
            composed = try
                product_basins(cib; structure = st)
            catch e
                e isa ArgumentError || rethrow()
                nothing
            end
            if composed === nothing
                result["composition_note"] =
                    "The islands' consistent scenarios combine to too many to " *
                    "enumerate in one table — analyse the islands separately " *
                    "(split_cib) and combine only the pieces you need."
            else
                total = Float64(composed.scenarioCount)
                compOrder = sortperm(composed.basinSizes; rev = true)
                result["composition"] = Dict(
                    "count" => length(composed.fixedPoints),
                    "cycles" => json_count(composed.cycleCount),
                    "scenarios" => [Dict(
                        "variants" => variant_names(cib, composed.fixedPoints[k]),
                        "basin_size" => json_count(composed.basinSizes[k]),
                        "basin_pct" => round(100 * Float64(composed.basinSizes[k]) / total;
                                             sigdigits = 4)) for k in compOrder])
            end
        else
            result["composition_note"] =
                "The model decomposes, but the largest island alone has " *
                "$(largestIsland) scenarios — analyse the islands separately " *
                "(estimate, or exact streaming via scripts/) and multiply."
        end
    end
    result["compute_s"] = round(compute; digits = 3)
    return result
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
            respond(sock, 400, "text/plain", "Run a basin analysis or estimate first.")
        else
            respond(sock, 200, "text/csv; charset=utf-8", state.last_csv;
                    extra = ["Content-Disposition" =>
                             "attachment; filename=\"" * state.last_csv_name * "\""])
        end
    elseif method == "POST" && path == "/analyze"
        mode = get(params, "mode", "consistent")
        try
            parse_s = @elapsed (cib = load_from_text(body))
            result = if mode == "basins"
                try
                    analyze_basins(cib, state)
                catch e
                    # find_basins refuses models whose tables cannot fit (an
                    # ArgumentError with directions). Surface that as guidance
                    # with the working alternatives, not as a failure.
                    e isa ArgumentError || rethrow()
                    Dict{String,Any}("mode" => "basins_too_big",
                                     "total" => json_count(scenario_count(cib)),
                                     "message" => e.msg,
                                     "threads" => Threads.nthreads())
                end
            elseif mode == "estimate"
                requested = something(tryparse(Int, get(params, "samples", "")), 1_000_000)
                analyze_estimate(cib, state; samples = clamp(requested, 1_000, 100_000_000))
            elseif mode == "structure"
                analyze_structure(cib)
            else
                analyze_consistent(cib)
            end
            result["parse_s"] = round(parse_s; digits = 3)
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
  select { font:inherit; border:1px solid var(--line); background:#fff;
           padding:8px 10px; border-radius:8px; color:var(--ink); }
  .guide { background:#fff8e6; border:1px solid #eadfa9; border-radius:8px;
           padding:12px 14px; white-space:pre-wrap; font-size:13px; }
  .cihint { color:var(--muted); font-size:12px; }
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
     scenarios exactly, analyse basins of attraction (exact where the space
     allows, estimated with error bars at any size), or map the model's
     influence structure.</p>
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
      <button class="green" id="btnBasins" disabled
              title="Walks every scenario — exact, needs the space to fit in memory">Exact Basins</button>
      <button class="green" id="btnEstimate" disabled
              title="Basin shares with confidence intervals — any size, seconds. The consistent scenarios themselves are still found exactly first">Estimate Basin Shares</button>
      <select id="samples" title="Sample count for the estimate">
        <option value="100000">100k samples</option>
        <option value="1000000" selected>1M samples</option>
        <option value="10000000">10M samples</option>
      </select>
      <button id="btnStructure" disabled
              title="Which descriptors actually influence which; independent islands and never-moving descriptors — the exact route for huge decomposable models">Structure</button>
      <span style="flex:1"></span>
      <button id="btnQuit" title="Stop the app">Quit</button>
    </div>
    <div id="status"></div>
  </div>
  <div class="card" id="results" style="display:none;">
    <div class="toolbar">
      <div class="summary" id="summary" style="flex:1"></div>
      <a id="exportLink" href="/export.csv" style="display:none;">
        <button class="green">Export CSV</button></a>
    </div>
    <div class="tablewrap"><table id="table"></table></div>
  </div>
  <footer>Cross-Impact Balances desktop app · runs locally on your machine</footer>
</main>
<script>
let scwText = null;
const $ = id => document.getElementById(id);
const RUN_BUTTONS = ["btnConsistent", "btnBasins", "btnEstimate", "btnStructure"];

$("file").addEventListener("change", e => {
  const f = e.target.files[0];
  if (!f) return;
  const reader = new FileReader();
  reader.onload = () => {
    scwText = reader.result;
    $("fname").textContent = f.name + " (" + f.size.toLocaleString() + " bytes)";
    RUN_BUTTONS.forEach(b => $(b).disabled = false);
    $("status").textContent = "";
  };
  reader.readAsText(f);
});

function setBusy(msg) {
  $("status").innerHTML = '<span class="spinner"></span>' + msg;
  RUN_BUTTONS.forEach(b => $(b).disabled = true);
}
function clearBusy() {
  $("status").textContent = "";
  RUN_BUTTONS.forEach(b => $(b).disabled = false);
}

const BUSY_TEXT = {
  consistent: "Searching for consistent scenarios…",
  basins: "Running exact basin analysis… (this can take a moment on large files)",
  estimate: "Estimating basin shares… (the consistent scenarios are found exactly first)",
  structure: "Mapping the influence structure…"
};

async function analyze(mode) {
  if (!scwText) return;
  setBusy(BUSY_TEXT[mode] || "Working…");
  try {
    const t0 = performance.now();
    let url = "/analyze?mode=" + mode;
    if (mode === "estimate") url += "&samples=" + $("samples").value;
    const res = await fetch(url, { method:"POST", body: scwText });
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

// Counts may arrive as strings (the server sends any count past 2^53 that
// way, because a JSON number would silently lose its last digits here in
// the browser). Strings pass through untouched — exact beats pretty.
function fmt(x) { return typeof x === "string" ? x : Number(x).toLocaleString(); }

function timingLine(data, secs) {
  let timing = `${secs}s`;
  if (data.compute_s !== undefined) {
    timing = `${data.compute_s}s compute on ${data.threads} thread`
           + (data.threads === 1 ? "" : "s")
           + (data.parse_s !== undefined ? ` · ${data.parse_s}s parse` : "")
           + ` · ${secs}s total`;
    if (data.threads === 1)
      timing += " — single-threaded! Launch via the Start-menu shortcut to use all cores";
  }
  return timing;
}

// One header row from column names; anything after `firstNumeric` is numeric.
function headerRow(cols, firstNumeric) {
  let html = "<thead><tr>";
  cols.forEach((c, i) => {
    const numeric = i >= firstNumeric || c === "Signature";
    html += `<th class="${numeric ? 'num' : ''}">${escapeHtml(c)}</th>`;
  });
  return html + "</tr></thead>";
}

function render(data, secs) {
  const timing = timingLine(data, secs);
  let sum = "";
  let html = "";
  let showExport = false;

  if (data.mode === "basins_too_big") {
    sum = `This model has <b>${fmt(data.total)}</b> scenarios — too large for the `
        + `exact table analysis on this machine. `
        + `<button class="green" id="btnEstimateNudge">Estimate basin shares instead</button>`
        + `<div class="guide" style="margin-top:10px;">${escapeHtml(data.message)}</div>`;
  } else if (data.mode === "structure") {
    const islands = data.islands.length;
    sum = `<b>${islands}</b> independent island${islands === 1 ? "" : "s"} across `
        + `${fmt(data.total)} scenarios · ${data.influences} active influence(s). `;
    if (data.dials.length)
      sum += `<br>Dials (never move; each setting can be analysed separately): `
           + data.dials.map(d => `<b>${escapeHtml(d.descriptor)}</b>`).join(", ") + ".";
    if (data.forced.length)
      sum += `<br>Jump-and-stay: ` + data.forced.map(d =>
             `<b>${escapeHtml(d.descriptor)}</b> settles into ${d.variants.map(escapeHtml).join("/")}`)
             .join("; ") + ".";
    if (data.composition) {
      sum += `<br>The islands compose: <b>${data.composition.count}</b> consistent `
           + `scenario(s) with <b>exact</b> basin sizes below; `
           + `${fmt(data.composition.cycles)} start(s) cycle.`;
      showExport = false;
      const comp = data.composition.scenarios;
      const maxBasin = Math.max(...comp.map(s => Number(String(s.basin_size).replace(/,/g,"")) || 0), 1);
      html = headerRow(["#", ...data.descriptors, "Basin size (exact)", "Basin %", ""],
                       1 + data.descriptors.length) + "<tbody>";
      comp.forEach((s, i) => {
        html += `<tr><td class="num">${i + 1}</td>`;
        s.variants.forEach(v => html += `<td>${escapeHtml(v)}</td>`);
        const w = Math.max(2, Math.round(70 * (Number(String(s.basin_size).replace(/,/g,"")) || 0) / maxBasin));
        html += `<td class="num">${fmt(s.basin_size)}</td><td class="num">${s.basin_pct}</td>`
              + `<td><div class="bar" style="width:${w}px"></div></td></tr>`;
      });
      html += "</tbody>";
    } else {
      if (data.composition_note)
        sum += `<div class="guide" style="margin-top:10px;">${escapeHtml(data.composition_note)}</div>`;
      html = headerRow(["#", "Island descriptors", "Scenarios"], 2) + "<tbody>";
      data.islands.forEach((isl, i) => {
        html += `<tr><td class="num">${i + 1}</td>`
              + `<td>${isl.descriptors.map(escapeHtml).join(", ")}</td>`
              + `<td class="num">${fmt(isl.scenarios)}</td></tr>`;
      });
      html += "</tbody>";
    }
  } else if (data.mode === "estimate") {
    sum = `<b>${data.count}</b> consistent scenario(s) — found <b>exactly</b> — with `
        + `shares estimated from ${fmt(data.samples)} samples over ${fmt(data.total)} `
        + `scenarios. Cycles: ${data.cycle_share_pct}% `
        + `<span class="cihint">[${data.cycle_ci_lo_pct}–${data.cycle_ci_hi_pct}]</span>. `
        + `Seed ${data.seed} — same seed, same result. `
        + `<span class="muted">(${timing})</span>`;
    showExport = true;
    const maxShare = Math.max(...data.scenarios.map(s => s.share_pct), 1e-9);
    html = headerRow(["#", ...data.descriptors, "Hits", "Share %", "95% CI", "Est. size", ""],
                     1 + data.descriptors.length) + "<tbody>";
    data.scenarios.forEach((s, i) => {
      html += `<tr><td class="num">${i + 1}</td>`;
      s.variants.forEach(v => html += `<td>${escapeHtml(v)}</td>`);
      if (s.hits === 0) {
        html += `<td class="num">0</td>`
              + `<td class="num muted">≤ ${s.ci_hi_pct}%</td>`
              + `<td class="num muted">upper bound</td><td class="num muted">—</td><td></td>`;
      } else {
        const w = Math.max(2, Math.round(70 * s.share_pct / maxShare));
        html += `<td class="num">${s.hits.toLocaleString()}</td>`
              + `<td class="num">${s.share_pct}</td>`
              + `<td class="num cihint">${s.ci_lo_pct}–${s.ci_hi_pct}</td>`
              + `<td class="num">${s.est_size.toLocaleString()}</td>`
              + `<td><div class="bar" style="width:${w}px"></div></td>`;
      }
      html += "</tr>";
    });
    html += "</tbody>";
  } else {
    const isBasins = data.mode === "basins";
    if (isBasins) {
      const cov = (100 * data.covered / data.total).toFixed(2);
      sum = `<b>${data.count}</b> consistent scenario(s) across a space of `
          + `${fmt(data.total)}. They attract <b>${cov}%</b> of all `
          + `scenarios; ${fmt(data.cycles)} start(s) fall into `
          + `non-fixed-point cycles. <span class="muted">(${timing})</span>`;
      showExport = true;
    } else {
      sum = `<b>${data.count}</b> consistent scenario(s) out of `
          + `${fmt(data.total)}. <span class="muted">(${timing})</span>`;
    }
    const maxBasin = isBasins && data.scenarios.length
        ? Math.max(...data.scenarios.map(s => s.basin_size)) : 0;
    const cols = ["#", "Signature", ...data.descriptors];
    if (isBasins) cols.push("Basin size", "Basin %", "");
    html = headerRow(cols, isBasins ? 2 + data.descriptors.length : cols.length) + "<tbody>";
    data.scenarios.forEach((s, i) => {
      html += `<tr><td class="num">${i + 1}</td>`
            + `<td class="num">${s.signature.toLocaleString()}</td>`;
      s.variants.forEach(v => html += `<td>${escapeHtml(v)}</td>`);
      if (isBasins) {
        html += `<td class="num">${s.basin_size.toLocaleString()}</td>`
              + `<td class="num">${s.basin_pct}</td>`;
        const w = maxBasin ? Math.max(2, Math.round(70 * s.basin_size / maxBasin)) : 0;
        html += `<td><div class="bar" style="width:${w}px"></div></td>`;
      }
      html += "</tr>";
    });
    html += "</tbody>";
  }

  $("summary").innerHTML = sum;
  $("table").innerHTML = html;
  $("exportLink").style.display = showExport ? "inline-block" : "none";
  $("results").style.display = "block";
  const nudge = $("btnEstimateNudge");
  if (nudge) nudge.addEventListener("click", () => analyze("estimate"));
}

$("btnConsistent").addEventListener("click", () => analyze("consistent"));
$("btnBasins").addEventListener("click", () => analyze("basins"));
$("btnEstimate").addEventListener("click", () => analyze("estimate"));
$("btnStructure").addEventListener("click", () => analyze("structure"));
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

#!/usr/bin/env julia
# Three-way CIB benchmark — JuCIB (Julia) leg.
#
# For each .scw file, times three modes on a PRE-PARSED model (find-only, so the
# parse cost is excluded and the numbers compare like-for-like with CIBSA):
#   1. full-enum, single-thread  (exhaustive=false, mc_threshold huge) — same
#      algorithm as CIBSA full enumeration; the apples-to-apples cell.
#   2. exhaustive, threaded      (exhaustive=true) — JuCIB's fast path / ground truth.
#   3. find_basins               — basin-of-attraction analysis (single-thread).
# Emits test/bench_results_julia.json.
#
# Run:  julia -t 10 --project=. test/bench_threeway.jl

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "sample_files")
const FILES  = ["CIB_global", "bench_medium", "bench_large",
                "bench_typical", "bench_xlarge", "bench_50x50"]

# warmup + n trials, return median seconds (n=2 -> faster of the two)
function timed(f; n::Int=3)
    f()                                   # warmup (JIT)
    ts = Float64[]
    for _ in 1:n
        t0 = time_ns(); f(); push!(ts, (time_ns() - t0) / 1e9)
    end
    sort!(ts)
    return ts[cld(n, 2)]
end

# parse only — passing a non-nothing kernel makes load_scw skip the search
preload(p) = load_scw(p; kernel = Vector{Vector{Int}}(), mc_threshold = 10^18)
canon(cib, k) = sort!([signature(cib, u) for u in k])

# ── minimal JSON writer (repo avoids JSON.jl) ───────────────────────────────
jstr(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
jval(x::Bool)           = x ? "true" : "false"
jval(x::Integer)        = string(x)
jval(x::AbstractFloat)  = (isnan(x) || isinf(x)) ? "null" : string(x)
jval(::Nothing)         = "null"
jval(x::AbstractString) = jstr(x)
jval(x::Vector{<:Pair}) = "{" * join([jstr(k) * ":" * jval(v) for (k, v) in x], ",") * "}"
jval(x::AbstractVector) = "[" * join([jval(v) for v in x], ",") * "]"

const NT = Threads.nthreads()
println("JuCIB three-way benchmark — Julia $(VERSION), threads = $NT")

records = Any[]

for name in FILES
    path = joinpath(SAMPLE, name * ".scw")
    if !isfile(path)
        println("  MISSING: $path — skipping"); continue
    end

    cib   = preload(path)
    total = max_signature(cib) + 1
    big   = total > 1_000_000
    ntr   = big ? 2 : 3
    do_fe1t = !big                         # single-thread full enum infeasible on huge spaces

    println("\n=== $name  ($total scenarios, $(cib.ndesc) descriptors) ===")

    parse_s = timed(() -> preload(path); n = ntr)

    fe1t_s = do_fe1t ? timed(() -> find_consistent(cib; exhaustive = false); n = ntr) : nothing

    exh_s = timed(() -> find_consistent(cib; exhaustive = true); n = ntr)
    kern  = find_consistent(cib; exhaustive = true)
    ksigs = canon(cib, kern)

    bas_s = timed(() -> find_basins(cib); n = ntr)
    fps, basins, cyc = find_basins(cib)
    fpsigs = [signature(cib, u) for u in fps]
    ord    = sortperm(fpsigs)
    basin_pairs = Any[ Pair{String,Any}["sig" => fpsigs[i], "size" => basins[i]] for i in ord ]
    min_basin   = isempty(basins) ? nothing : minimum(basins)

    verified = all(succession_step(cib, u) == u for u in kern)

    push!(records, Pair{String,Any}[
        "file" => name, "tool" => "julia", "threads" => NT,
        "total_scenarios" => total,
        "parse_s"         => round(parse_s, digits = 6),
        "full_enum_1t_s"  => fe1t_s === nothing ? nothing : round(fe1t_s, digits = 6),
        "exhaustive_s"    => round(exh_s, digits = 6),
        "basins_s"        => round(bas_s, digits = 6),
        "kernel_size"     => length(kern),
        "kernel_sigs"     => ksigs,
        "cycle_count"     => cyc,
        "min_basin"       => min_basin,
        "basins"          => basin_pairs,
        "verified_fixed_points" => verified,
    ])

    fe_str = fe1t_s === nothing ? "—(skip)" : "$(round(fe1t_s, digits = 4))s"
    println("  parse=$(round(parse_s, digits = 4))s  fe1t=$fe_str  " *
            "exh=$(round(exh_s, digits = 4))s  basins=$(round(bas_s, digits = 4))s")
    println("  kernel=$(length(kern))  cycles=$cyc  min_basin=$min_basin  verified=$verified")
    println("  sigs=$ksigs")
end

open(joinpath(@__DIR__, "bench_results_julia.json"), "w") do io
    write(io, jval(records))
end
println("\nWrote test/bench_results_julia.json  ($(length(records)) files)")

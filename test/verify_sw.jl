#!/usr/bin/env julia
# Verify ScenarioWizard's exported .sl solutions against JuCIB.
# Loads each <name>_sw.sl via load_scw(scw; sl_file=...), reduces SW's kernel to
# sorted signatures (JuCIB's space), and checks every SW scenario is a true fixed
# point (succession_step(u) == u). Emits test/bench_results_sw.json.
#
# Run:  julia --project=. test/verify_sw.jl

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "sample_files")

# Observed GUI wall-clock (coarse; screenshot-paced). The five small/mid files
# returned with no progress dialog (instant); only 60M showed a "Computing…" bar.
const SW_TIME = Dict(
    "CIB_global"   => (nothing, "instant (<~0.5 s, no progress dialog)"),
    "bench_medium" => (nothing, "instant (<~0.5 s, no progress dialog)"),
    "bench_large"  => (nothing, "instant (<~0.5 s, no progress dialog)"),
    "bench_typical"=> (nothing, "instant (<~0.5 s, no progress dialog)"),
    "bench_xlarge" => (nothing, "instant (<~0.5 s, no progress dialog)"),
    "bench_50x50"  => (5.0,     "≈2–6 s (progress dialog observed)"),
)

jstr(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
jval(x::Bool)           = x ? "true" : "false"
jval(x::Integer)        = string(x)
jval(x::AbstractFloat)  = (isnan(x) || isinf(x)) ? "null" : string(x)
jval(::Nothing)         = "null"
jval(x::AbstractString) = jstr(x)
jval(x::Vector{<:Pair}) = "{" * join([jstr(k) * ":" * jval(v) for (k, v) in x], ",") * "}"
jval(x::AbstractVector) = "[" * join([jval(v) for v in x], ",") * "]"

files = ["CIB_global", "bench_medium", "bench_large",
         "bench_typical", "bench_xlarge", "bench_50x50"]

records = Any[]
println("ScenarioWizard solution verification (via JuCIB):\n")

for name in files
    scw = joinpath(SAMPLE, name * ".scw")
    sl  = joinpath(SAMPLE, name * "_sw.sl")
    if !isfile(sl)
        println("  $name: no _sw.sl — skipping"); continue
    end
    cib = load_scw(scw; sl_file = sl)            # cib.kernel = SW solutions (0-based)
    sigs = sort([signature(cib, u) for u in cib.kernel])
    verified = all(succession_step(cib, u) == u for u in cib.kernel)
    t, note = SW_TIME[name]

    push!(records, Pair{String,Any}[
        "file" => name, "tool" => "scenariowizard",
        "kernel_size" => length(cib.kernel),
        "kernel_sigs" => sigs,
        "verified_fixed_points" => verified,
        "time_s_approx" => t,
        "time_note" => note,
    ])
    flag = verified ? "OK" : "*** NOT ALL FIXED POINTS ***"
    println("  $name: count=$(length(cib.kernel))  verified=$verified [$flag]")
    println("    sigs=$sigs")
end

open(joinpath(@__DIR__, "bench_results_sw.json"), "w") do io
    write(io, jval(records))
end
println("\nWrote test/bench_results_sw.json  ($(length(records)) files)")

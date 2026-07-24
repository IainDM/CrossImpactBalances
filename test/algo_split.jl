#!/usr/bin/env julia
# Separate the algorithmic speedup from the threading speedup by timing the
# exhaustive fixed-point sweep at the CURRENT thread count.
# Run once with -t 1 and once with -t 10; the ratio is the threading speedup.
using CrossImpactBalances
const SAMPLE = joinpath(@__DIR__, "sample_files")
timed(f; n=3) = (f(); ts=Float64[]; for _ in 1:n; t=time_ns(); f(); push!(ts,(time_ns()-t)/1e9); end; sort!(ts); ts[cld(n,2)])
preload(p) = load_scw(p; kernel=Vector{Vector{Int}}())
println("threads = ", Threads.nthreads())
for name in ["bench_typical", "bench_xlarge", "bench_50x50"]
    cib = preload(joinpath(SAMPLE, name * ".scw"))
    big = (max_signature(cib) + 1) > 1_000_000
    swp = timed(() -> find_consistent(cib; algorithm=:sweep); n = big ? 2 : 3)
    bnb = timed(() -> find_consistent(cib; algorithm=:bnb);   n = big ? 2 : 3)
    println("$name: sweep=$(round(swp, digits=6))s  bnb=$(round(bnb, digits=6))s")
end

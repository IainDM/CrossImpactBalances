#!/usr/bin/env julia
# Separate the algorithmic speedup from the threading speedup by timing the
# exhaustive fixed-point sweep at the CURRENT thread count. Compare against the
# single-thread succession-walk (find_consistent exhaustive=false) which is the
# "same algorithm as CIBSA".  Run once with -t 1 and once with -t 10.
using CrossImpactBalances
const SAMPLE = joinpath(@__DIR__, "sample_files")
timed(f; n=3) = (f(); ts=Float64[]; for _ in 1:n; t=time_ns(); f(); push!(ts,(time_ns()-t)/1e9); end; sort!(ts); ts[cld(n,2)])
preload(p) = load_scw(p; kernel=Vector{Vector{Int}}(), mc_threshold=10^18)
println("threads = ", Threads.nthreads())
for name in ["bench_typical", "bench_xlarge", "bench_50x50"]
    cib = preload(joinpath(SAMPLE, name * ".scw"))
    big = (max_signature(cib) + 1) > 1_000_000
    exh  = timed(() -> find_consistent(cib; exhaustive=true);  n = big ? 2 : 3)
    walk = big ? nothing : timed(() -> find_consistent(cib; exhaustive=false); n = 3)
    println("$name: exhaustive=$(round(exh, digits=6))s  walk_1t=$(walk===nothing ? "skip" : string(round(walk, digits=6))*"s")")
end

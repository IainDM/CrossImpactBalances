#!/usr/bin/env julia
# Clean average of find_basins wall-clock on bench_50x50 (60M scenarios).
# Reads the .scw from a fixed absolute path so both branches use identical input.
using CrossImpactBalances
const SAMPLE = raw"D:\GitHub\CrossImpactBalances\test\sample_files"
const N = 3
cib = load_scw(joinpath(SAMPLE, "bench_50x50.scw"); kernel=Vector{Vector{Int}}())
find_basins(cib)  # warmup (exclude JIT)
ts = Float64[]
for i in 1:N
    t = time_ns(); find_basins(cib); dt = (time_ns() - t) / 1e9
    push!(ts, dt); println("  trial $i: $(round(dt, digits=2)) s")
end
m = sum(ts) / length(ts)
println("MEAN = $(round(m, digits=2)) s   (min $(round(minimum(ts),digits=2)), max $(round(maximum(ts),digits=2)), n=$N, threads=$(Threads.nthreads()))")
println("pkg  = ", pkgdir(CrossImpactBalances))

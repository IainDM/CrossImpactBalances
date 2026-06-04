"""
Threaded exhaustive search on `bench_typical.scw` (10 descriptors × 3 variants
= 59,049 scenarios). Demonstrates the `exhaustive=true` path.

Start Julia with multiple threads to see the speedup:
    julia --project=. -t auto examples/03_threaded_exhaustive.jl
"""

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files", "bench_typical.scw")

println("Threads in use: ", Threads.nthreads())
println()

t0 = time_ns()
cib = load_scw(SAMPLE; exhaustive=true)
elapsed = (time_ns() - t0) / 1e9
println("Exhaustive find_consistent: $(round(elapsed, digits=3))s")
println("Fixed points found: $(length(cib.kernel))")
for u in cib.kernel
    println("  sig=$(signature(cib, u))")
end
println()

t0 = time_ns()
fps, basins, cyc = find_basins(cib)
elapsed = (time_ns() - t0) / 1e9
println("Basin analysis: $(round(elapsed, digits=3))s")
for (i, u) in enumerate(fps)
    println("  fp sig=$(signature(cib, u))  basin=$(basins[i])")
end
println("  scenarios in cycles: $cyc")

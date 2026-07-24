#!/usr/bin/env julia
# Compare two checkouts of CrossImpactBalances on identical input files.
# Reads .scw from a FIXED absolute path so both branches benchmark the same bytes.
# Run from each worktree with the same `-t`:
#   julia -t 10 --project=. <thisfile>
using CrossImpactBalances
const SAMPLE = raw"D:\GitHub\CrossImpactBalances\test\sample_files"
timed(f; n=3) = (f(); ts=Float64[]; for _ in 1:n; t=time_ns(); f(); push!(ts,(time_ns()-t)/1e9); end; sort!(ts); ts[cld(n,2)])
preload(p) = load_scw(p; kernel=Vector{Vector{Int}}())
println("pkg     = ", pkgdir(CrossImpactBalances))
println("threads = ", Threads.nthreads())
for name in ["bench_typical", "bench_xlarge", "bench_50x50"]
    cib = preload(joinpath(SAMPLE, name * ".scw"))
    n   = max_signature(cib) + 1
    ntr = n > 1_000_000 ? 2 : 3
    exh  = timed(() -> find_consistent(cib); n=ntr)
    bas  = timed(() -> find_basins(cib); n=ntr)
    k    = length(find_consistent(cib))
    println("$name (n=$n, kernel=$k):  find_consistent=$(round(exh,digits=6))s  basins=$(round(bas,digits=6))s")
end

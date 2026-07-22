# Determinism / race check for the threaded `find_basins`: run it many times and assert
# every run returns the IDENTICAL (fixed points, basin sizes, cycle count). Prints a
# fingerprint so runs at different `-t` can be compared — they MUST match, since the
# result is defined to be independent of thread count.
#
#   julia --project=. -t 8 bench/determinism.jl [nested|50x50] [runs]
#   # then compare the fingerprint against:  julia --project=. -t 1 bench/determinism.jl nested 1
using CrossImpactBalances
const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files")
preload(p) = load_scw(p; kernel = Vector{Vector{Int}}())

model = get(ARGS, 1, "nested")
runs = parse(Int, get(ARGS, 2, "8"))
cib = preload(joinpath(SAMPLE, model == "nested" ? "CIB_nested.scw" : "bench_50x50.scw"))
n = max_signature(cib) + 1

ref = nothing
for i in 1:runs
    fps, sizes, cyc = find_basins(cib)
    key = ([signature(cib, u) for u in fps], sizes, cyc)   # canonical, order-independent
    @assert sum(sizes) + cyc == n "run $i: coverage broken"
    if ref === nothing
        global ref = key
    else
        @assert key == ref "run $i DIFFERS from run 1 — RACE / nondeterminism!"
    end
end
sigs, sizes, cyc = ref
println(model, " -t", Threads.nthreads(), " x", runs, " runs ALL IDENTICAL")
println("  nfp=", length(sigs), " cyc=", cyc, " sum(basins)=", sum(sizes), " fingerprint=", hash(ref))

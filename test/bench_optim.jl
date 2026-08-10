# Post-optimization benchmark: exhaustive sweep vs branch-and-bound vs basin
# analysis on the shipped sample files, with a cross-algorithm agreement check.
#
#   julia -t auto --project=. test/bench_optim.jl
#
# Times are find-only medians of 3 (a warm-up run is discarded); parsing is
# excluded. B&B nodes = partial scenarios expanded by the pruned search,
# reported against the full space size.

using CrossImpactBalances

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")

med3(f) = (f(); sort([(@elapsed f()) for _ in 1:3])[2])

println("Julia $(VERSION), threads = $(Threads.nthreads())")
println()
println("| File | Scenarios | sweep (s) | B&B (s) | B&B nodes visited | find_basins (s) |")
println("|---|---:|---:|---:|---:|---:|")

for name in ["bench_typical", "bench_xlarge", "bench_50x50"]
    cib = load_scw(joinpath(SAMPLE_DIR, "$name.scw"); kernel=Vector{Vector{Int}}())
    n = max_signature(cib) + 1

    t_sweep = med3(() -> find_consistent(cib; algorithm=:sweep))
    t_bnb   = med3(() -> find_consistent(cib; algorithm=:bnb))
    t_bas   = med3(() -> find_basins(cib))

    # _bnb_search is the whole branch-and-bound path, descriptor ordering
    # included — the same thing find_consistent(:bnb) above just timed.
    _, nodes = CrossImpactBalances._bnb_search(cib; node_budget=typemax(Int))

    ks = [signature(cib, u) for u in find_consistent(cib; algorithm=:sweep)]
    kb = [signature(cib, u) for u in find_consistent(cib; algorithm=:bnb)]
    fps, _, _ = find_basins(cib)
    kf = [signature(cib, u) for u in fps]
    ks == kb == kf || error("algorithms disagree on $name: sweep=$ks bnb=$kb basins=$kf")

    pct = round(100 * nodes / n; sigdigits=3)
    println("| $name | $n | $(round(t_sweep; sigdigits=3)) | ",
            "$(round(t_bnb; sigdigits=3)) | $nodes ($pct%) | ",
            "$(round(t_bas; sigdigits=3)) |")
end

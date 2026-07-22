# Per-stage wall-time breakdown of `find_basins`: the parallel successor-table build
# (`_successor_table!`) vs the resolve+tally (`_resolve_and_tally`). Useful for seeing
# where time goes and how each stage scales with threads.
#
#   julia --project=. -t 8 bench/stage_breakdown.jl [nested|50x50|typical]
#   for t in 1 2 4 8; do julia --project=. -t $t bench/stage_breakdown.jl nested; done
using CrossImpactBalances
const CIB = CrossImpactBalances
const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files")
preload(p) = load_scw(p; kernel = Vector{Vector{Int}}())

const FILES = Dict("nested" => "CIB_nested.scw", "50x50" => "bench_50x50.scw",
                   "typical" => "bench_typical.scw")

function breakdown(cib; label = "", trials = 3)
    n = max_signature(cib) + 1
    T = CIB._score_type(cib)
    S = n <= Int(typemax(Int32)) - 1 ? Int32 : Int64
    cimT = Matrix{T}(cib.cim_t)
    best = (Inf, Inf)
    local fps, sizes, cyc
    for _ in 1:trials
        succ = Vector{S}(undef, n)
        GC.gc(); t0 = time_ns(); CIB._successor_table!(succ, cib, cimT); tb = (time_ns() - t0) / 1e9
        GC.gc(); t0 = time_ns(); fps, sizes, cyc = CIB._resolve_and_tally(succ, n); tr = (time_ns() - t0) / 1e9
        tb + tr < sum(best) && (best = (tb, tr))
    end
    tb, tr = best; tot = tb + tr
    println(rpad(label, 8), " n=", rpad(n, 10), " threads=", Threads.nthreads(),
            " | table=", rpad(round(tb, digits=3), 7), "s(", lpad(round(Int, 100tb/tot), 3), "%)",
            " resolve+tally=", rpad(round(tr, digits=3), 7), "s(", lpad(round(Int, 100tr/tot), 3), "%)",
            " | wall=", round(tot, digits=3), "s")
    return fps, sizes, cyc
end

# Warm/compile on the small model, then break down the requested one.
breakdown(preload(joinpath(SAMPLE, FILES["typical"])); label = "warmup", trials = 1)
model = get(ARGS, 1, "50x50")
cib = preload(joinpath(SAMPLE, get(FILES, model, FILES["50x50"])))
fps, sizes, cyc = breakdown(cib; label = model)
@assert sum(sizes) + cyc == max_signature(cib) + 1 "coverage check failed"

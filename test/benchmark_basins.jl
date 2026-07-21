"""
Basin-analysis benchmark: times `find_basins` on the large sample models and
verifies the results against pinned references.

Run with:
    julia --project=. -t 4 test/benchmark_basins.jl            # bench_50x50 only
    julia --project=. -t 4 test/benchmark_basins.jl nested     # + CIB_nested (408M states)

Timings are min/median of `--trials` runs (default 3) after a warmup on
bench_typical (which compiles the same method instances the large models use).
"""

using CrossImpactBalances

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")

# Skip the kernel search at load time; find_basins recomputes fixed points.
preload(path) = load_scw(path; kernel=Vector{Vector{Int}}())

# Pinned references (test/bench_results_julia.json, verified against Python/SW).
const EXPECT_50X50 = Dict(1984298 => 3, 8693546 => 578, 11540522 => 113420,
                          19807290 => 95, 31030843 => 5633)
const CYC_50X50 = 60_346_447
const EXPECT_TYPICAL = Dict(13785 => 52329, 13839 => 1539)
const CYC_TYPICAL = 5181

function bench(name, path; expect=nothing, cyc_expect=nothing, trials=3)
    cib = preload(path)
    n = max_signature(cib) + 1
    ts = Float64[]
    local fps, basins, cyc
    for _ in 1:trials
        GC.gc()
        t0 = time_ns()
        fps, basins, cyc = find_basins(cib)
        push!(ts, (time_ns() - t0) / 1e9)
    end
    sort!(ts)
    sigs = [signature(cib, u) for u in fps]
    @assert sum(basins) + cyc == n "$name: coverage check failed (sum+cyc=$(sum(basins)+cyc) != $n)"
    @assert all(succession_step(cib, u) == u for u in fps) "$name: non-fixed-point returned"
    if expect !== nothing
        @assert Dict(zip(sigs, basins)) == expect "$name: fp sigs/basin sizes mismatch"
        @assert cyc == cyc_expect "$name: cycle_count $cyc != expected $cyc_expect"
    end
    med = ts[cld(length(ts), 2)]
    println(rpad(name, 24), " n=", rpad(n, 11), " min=", rpad("$(round(ts[1], digits=3))s", 9),
            " median=", rpad("$(round(med, digits=3))s", 9),
            " nfp=", length(fps), " cyc=", cyc,
            " maxrss=", round(Sys.maxrss() / 2^30, digits=2), "GiB")
    return ts[1]
end

println("Julia $(VERSION), $(Threads.nthreads()) threads")

# Warmup + small-model sanity (also compiles the Int16 fast path)
bench("bench_typical (warmup)", joinpath(SAMPLE_DIR, "bench_typical.scw");
      expect=EXPECT_TYPICAL, cyc_expect=CYC_TYPICAL, trials=1)

bench("bench_50x50", joinpath(SAMPLE_DIR, "bench_50x50.scw");
      expect=EXPECT_50X50, cyc_expect=CYC_50X50)

if "nested" in ARGS
    path = joinpath(SAMPLE_DIR, "CIB_nested.scw")
    cib = preload(path)
    sl_sigs = Set(signature(cib, u) for u in
                  load_solutions(cib, joinpath(SAMPLE_DIR, "CIB_nested.sl")))
    t0 = time_ns()
    fps, basins, cyc = find_basins(cib)
    t = (time_ns() - t0) / 1e9
    n = max_signature(cib) + 1
    @assert sum(basins) + cyc == n "CIB_nested: coverage check failed"
    @assert Set(signature(cib, u) for u in fps) == sl_sigs "CIB_nested: fps != .sl solutions"
    println(rpad("CIB_nested", 24), " n=", rpad(n, 11), " time=", round(t, digits=3),
            "s  nfp=", length(fps), " cyc=", cyc,
            " maxrss=", round(Sys.maxrss() / 2^30, digits=2), "GiB")
end

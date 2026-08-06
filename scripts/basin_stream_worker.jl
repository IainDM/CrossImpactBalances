# One machine's share of an exact streaming basin analysis.
#
#   julia --project=. -t auto scripts/basin_stream_worker.jl model.scw FIRST LAST out.tsv [cache_bytes]
#
# Runs find_basins(...; method=:stream, signature_range=FIRST:LAST) — exact,
# memory-flat — and writes the partial counts as a small TSV. Split the space
# 0:(scenario_count-1) into disjoint ranges, run one worker per range (SLURM
# array jobs, a lab of desktops, anything), then combine the TSVs with
# scripts/basin_stream_merge.jl. Ranges may be any sizes as long as together
# they tile the space exactly; the merge refuses gaps and overlaps.
#
# The TSV is line-oriented and self-describing:
#   M <model file name>            which model this partial belongs to
#   R <first> <last>               the start range this worker covered
#   F <fixed point signature> <count>     one line per fixed point seen
#   C <cycle count>                starts in this range that end in a cycle
#
# Fixed points are recorded by SIGNATURE (Int64 — the guard in max_signature
# has already ensured signatures fit, or this analysis could not run at all).
#
# Measure before you commit a cluster: bench/stream_calibration.jl reports
# this machine's starts/second and extrapolates the wall-clock for a target
# space. At ~10^12–10^13 starts, exact streaming is CPU-days spread across
# your cores; beyond that, use estimate_basins / product_basins instead.

using CrossImpactBalances

function main(args)
    if length(args) < 4
        println(stderr, "usage: basin_stream_worker.jl model.scw FIRST LAST out.tsv [cache_bytes]")
        return 2
    end
    scwPath = args[1]
    firstSignature = parse(Int, args[2])
    lastSignature = parse(Int, args[3])
    outPath = args[4]
    cacheBytes = length(args) >= 5 ? parse(Int, args[5]) : 1 << 30

    cib = load_scw(scwPath; kernel=Vector{Vector{Int}}())   # parse only; no kernel search
    println(stderr, "basin_stream_worker: ", basename(scwPath), " — ",
            scenario_count(cib), " scenarios total; this worker takes ",
            firstSignature, ":", lastSignature, " on ", Threads.nthreads(), " threads")

    startedAt = time()
    fixedPoints, basinSizes, cycleCount =
        find_basins(cib; method=:stream, signature_range=firstSignature:lastSignature,
                    cache_bytes=cacheBytes, progress=true)
    elapsed = time() - startedAt

    # Write atomically: a half-written partial must never look mergeable.
    temporary = outPath * ".tmp"
    open(temporary, "w") do io
        println(io, "M\t", basename(scwPath))
        println(io, "R\t", firstSignature, "\t", lastSignature)
        for (scenario, size) in zip(fixedPoints, basinSizes)
            println(io, "F\t", signature(cib, scenario), "\t", size)
        end
        println(io, "C\t", cycleCount)
    end
    mv(temporary, outPath; force=true)

    covered = sum(basinSizes; init=0) + cycleCount
    expected = lastSignature - firstSignature + 1
    covered == expected || error("coverage check failed: $covered != $expected")
    println(stderr, "basin_stream_worker: done in ", round(elapsed; digits=1), "s — ",
            length(fixedPoints), " fixed points seen, ", cycleCount,
            " cycle starts; wrote ", outPath)
    return 0
end

exit(main(ARGS))

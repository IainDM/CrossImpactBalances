# Measure this machine's streaming-basins throughput and forecast wall-clock
# for a target scenario space BEFORE committing cluster time.
#
#   julia --project=. -t 8 bench/stream_calibration.jl [nested|50x50|typical] [starts] [full]
#
# Streams `starts` starting scenarios (default 10^7) of the chosen model with
# the memo cache on and off, reports starts/second, and extrapolates the
# wall-clock for the two ZIRIUS-scale spaces (6.7e12 and 3.55e18 scenarios)
# at several core counts, assuming linear scaling from this machine's
# per-core rate — real clusters come close for this workload, since workers
# share nothing but a counter.
#
# `full` additionally streams the model's ENTIRE space and asserts exact
# equality against the table method — the end-to-end proof on a real model
# (408M scenarios for `nested`: expect a few GiB of RAM for the table side
# and a stream measured in minutes, not seconds).

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files")
const FILES = Dict("nested" => "CIB_nested.scw", "50x50" => "bench_50x50.scw",
                   "typical" => "bench_typical.scw")
preload(p) = load_scw(p; kernel = Vector{Vector{Int}}())

function measure(cib, starts; cache_bytes)
    range = 0:min(starts, max_signature(cib) + 1) - 1
    t0 = time_ns()
    fps, sizes, cyc = find_basins(cib; method=:stream, signature_range=range,
                                  cache_bytes=cache_bytes)
    elapsed = (time_ns() - t0) / 1e9
    @assert sum(sizes; init=0) + cyc == length(range) "coverage check failed"
    return elapsed, length(range) / elapsed
end

function forecast(rate_per_core, label, n)
    print(rpad(label, 22), " ")
    for cores in (Threads.nthreads(), 32, 128, 512, 2048)
        seconds = Float64(n) / (rate_per_core * cores)
        text = seconds < 3600 ? string(round(seconds / 60, digits=1), "m") :
               seconds < 86400 ? string(round(seconds / 3600, digits=1), "h") :
               seconds < 86400 * 365 * 3 ? string(round(seconds / 86400, digits=1), "d") :
               string(round(seconds / (86400 * 365), sigdigits=2), "y")
        print(lpad(string(cores, " cores: ", text), 19))
    end
    println()
end

model = get(ARGS, 1, "nested")
starts = length(ARGS) >= 2 && tryparse(Int, ARGS[2]) !== nothing ? parse(Int, ARGS[2]) : 10^7

# Warmup compiles the same method instances the measured model uses.
let warm = preload(joinpath(SAMPLE, FILES["typical"]))
    find_basins(warm; method=:stream)
end

cib = preload(joinpath(SAMPLE, get(FILES, model, FILES["nested"])))
println("Julia $(VERSION), $(Threads.nthreads()) threads — ", model, ", ",
        scenario_count(cib), " scenarios, measuring ", starts, " starts")

tCached, rateCached = measure(cib, starts; cache_bytes=1 << 30)
tPlain, ratePlain = measure(cib, starts; cache_bytes=0)
println("with 1 GiB memo cache : ", round(tCached, digits=2), "s  = ",
        round(Int, rateCached), " starts/s (", round(Int, rateCached / Threads.nthreads()),
        " per core)")
println("without cache         : ", round(tPlain, digits=2), "s  = ",
        round(Int, ratePlain), " starts/s (", round(Int, ratePlain / Threads.nthreads()),
        " per core)")

perCore = rateCached / Threads.nthreads()
println("\nWall-clock forecast at this per-core rate (linear scaling):")
forecast(perCore, "6.7e12 (2^22*3^13)", Int128(6_687_075_336_192))
forecast(perCore, "3.55e18 (2^22*3^25)", Int128(2)^22 * Int128(3)^25)
println("\nThe 3.55e18 row is the point: exact enumeration stays out of reach at any")
println("plausible core count — that space needs estimate_basins (shares with error")
println("bars) and influence_structure/product_basins (exact, if the matrix splits).")

if "full" in ARGS
    println("\nfull: streaming the whole ", max_signature(cib) + 1, "-scenario space...")
    t0 = time_ns()
    streamed = find_basins(cib; method=:stream, progress=true)
    tStream = (time_ns() - t0) / 1e9
    t0 = time_ns()
    tabled = find_basins(cib; method=:table)
    tTable = (time_ns() - t0) / 1e9
    @assert streamed == tabled "STREAM AND TABLE DISAGREE — do not trust either"
    println("full: stream == table exactly (", length(tabled[1]), " fixed points, ",
            tabled[3], " cycle starts) — stream ", round(tStream, digits=1),
            "s vs table ", round(tTable, digits=1), "s")
end

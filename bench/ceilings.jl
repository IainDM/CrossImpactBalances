# Machine memory/compute ceilings, for interpreting `find_basins` numbers on a given box:
#   * STREAM triad  -> sustained DRAM bandwidth (GB/s)
#   * random chase a[a[i]] -> dependent-load latency (ns/access), the resolve's primitive
# Compare the resolve's ns/state against the chase latency to tell "cache-served" from
# "DRAM-latency-bound", and the table build's write bandwidth against STREAM.
#
#   julia --project=. -t 8 bench/ceilings.jl
using Base.Threads

function triad!(a, b, c, s)
    @threads :static for i in eachindex(a)
        @inbounds a[i] = b[i] + s * c[i]
    end
end
function stream_bw(gb)
    N = round(Int, gb * 2^30 / 8)
    a = fill(1.0, N); b = fill(2.0, N); c = fill(3.0, N)
    triad!(a, b, c, 1.5)
    best = Inf
    for _ in 1:5
        t0 = time_ns(); triad!(a, b, c, 1.5); best = min(best, (time_ns() - t0) / 1e9)
    end
    return 3 * 8 * N / best / 1e9      # 2 read + 1 write (ignores write-allocate)
end
function chase(gb, steps)
    N = round(Int, gb * 2^30 / 4)
    a = Vector{Int32}(undef, N)
    @inbounds for i in 1:N; a[i] = Int32((i * 2654435761) % N); end
    cur = Int32(0)
    a[1]                               # touch
    t0 = time_ns()
    @inbounds for _ in 1:steps; cur = a[cur + 1]; end
    dt = (time_ns() - t0) / 1e9
    return dt / steps * 1e9, cur       # ns/access
end

println("threads=", Threads.nthreads())
println("STREAM triad (2 GB) : ", round(stream_bw(2.0), digits=1), " GB/s")
ns, _ = chase(1.6, 50_000_000)
println("random chase (1.6 GB, 1 thread): ", round(ns, digits=1), " ns/access")

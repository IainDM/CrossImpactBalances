# ── Very large scenario spaces ──────────────────────────────────────────────
#
# Real CIB matrices reach scenario spaces of 10^12 and beyond — sizes where
# find_basins' tables cannot exist (they would need terabytes) and where even
# counting the space needs care (Int64 wraps past ~9.2e18; Float64 display
# mangles counts past ~9e15). This example walks the toolkit for that regime:
#
#   scenario_count       exact Int128 size, at any scale
#   find_basins          refuses politely when the tables cannot fit...
#     method=:stream     ...or computes exactly with flat memory, restrictable
#                        to a signature_range for splitting across machines
#   estimate_basins      basin SHARES with confidence intervals, in seconds,
#                        at any scale — the kernel itself stays exact
#   influence_structure  who actually influences whom; independent islands
#   product_basins       exact basins for a huge space, composed from islands
#
# Everything here runs in seconds on synthetic models; swap in your own .scw.

using CrossImpactBalances
using Random

# A model SHAPED like a real monster: 11 four-variant + 13 three-variant
# descriptors = 4^11 · 3^13 = 6,687,075,336,192 scenarios. (Real matrices of
# this size exist; this one is random, standing in for the shape only.)
rng = MersenneTwister(2026)
nvariants = vcat(fill(4, 11), fill(3, 13))
ndim = sum(nvariants)
cim = rand(rng, -3:3, ndim, ndim)
offset = 0
for count in nvariants                      # standard form: a descriptor casts
    cim[offset+1:offset+count, offset+1:offset+count] .= 0   # no votes on itself
    global offset += count
end
descriptors = ["D$i" for i in 1:length(nvariants)]
variants = Dict(descriptors[i] => ["V$(i)_$j" for j in 1:nvariants[i]]
                for i in eachindex(nvariants))
monster = CIB(descriptors, variants, nvariants, cim, permutedims(cim), ndim,
              length(nvariants), Vector{Vector{Int}}(),
              cumsum(vcat(0, nvariants[1:end-1])))

println("scenario_count: ", scenario_count(monster), " scenarios")

# find_basins knows it cannot build tables for this and says what to do
# instead — it never crashes, and never silently starts a multi-day job.
try
    find_basins(monster)
catch declined
    println("\nfind_basins declined, with directions:\n", declined.msg)
end

# Estimation works immediately, at any scale: uniform starts, exact walks,
# shares with Wilson intervals. The kernel would normally come from
# find_consistent (exact even here — branch-and-bound never enumerates);
# for this random stand-in we skip kernel discovery and sample cold.
estimate = estimate_basins(monster; samples=2_000, kernel=Vector{Vector{Int}}())
println("\nestimate_basins on the monster (2,000 samples):")
show(stdout, MIME"text/plain"(), estimate)
println()

# Exact streaming is the third option: flat memory, every start walked. On a
# space this size that is a cluster job (see bench/stream_calibration.jl and
# scripts/basin_stream_worker.jl); a small signature_range shows the shape.
fps, sizes, cycles = find_basins(monster; method=:stream, signature_range=0:99_999)
println("\nstream over starts 0:99999 — ", length(fps), " fixed points reached, ",
        cycles, " cycle starts, coverage ", sum(sizes; init=0) + cycles, "/100000")

# ── Decomposition: the exact route ──────────────────────────────────────────
# Two independent sub-models written into one matrix (block-diagonal
# judgments), plus one descriptor nobody votes on — a dial.
nv2 = [3, 3, 2, 2, 3]
dim2 = sum(nv2)
off2 = cumsum(vcat(0, nv2[1:end-1]))
cim2 = zeros(Int, dim2, dim2)
cim2[1:6, 1:6] .= rand(rng, -3:3, 6, 6)          # island one: D1, D2
cim2[7:10, 7:10] .= rand(rng, -3:3, 4, 4)        # island two: D3, D4
offset = 0                                        # D5 (columns 11:13): silent
for count in nv2
    cim2[offset+1:offset+count, offset+1:offset+count] .= 0
    global offset += count
end
d2 = ["D$i" for i in 1:5]
small = CIB(d2, Dict(d2[i] => ["V$(i)_$j" for j in 1:nv2[i]] for i in 1:5),
            nv2, cim2, permutedims(cim2), dim2, 5, Vector{Vector{Int}}(), off2)

structure = influence_structure(small)
println("\nThe influence map sees the seams:")
show(stdout, MIME"text/plain"(), structure)
println()

# product_basins runs find_basins per island and composes exactly — basin
# sizes multiply across islands, cycles compose, everything stays exact. On a
# genuinely huge decomposable model this is the only exact analysis there is.
composed = product_basins(small; structure=structure)
println("\nExact basins by composition:")
show(stdout, MIME"text/plain"(), composed)
println()

# A dial can also be pinned by hand — conditional analysis of one slice:
pinned = fix_descriptor(small, "D5", 0)
slice_fps, slice_sizes, slice_cycles = find_basins(pinned)
println("\nWith D5 pinned to its first variant: ", length(slice_fps),
        " consistent scenarios over ", scenario_count(pinned), " scenarios")

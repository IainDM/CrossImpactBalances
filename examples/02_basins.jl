"""
Basin-of-attraction analysis on `CIB_global`. Shows how many of the 36
starting scenarios converge to each consistent scenario, and how many fall
into non-fixed-point cycles.

Run from the repo root:
    julia --project=. examples/02_basins.jl
"""

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files", "CIB_global.scw")

cib = load_scw(SAMPLE)
fixed_points, basins, cycle_count = find_basins(cib)

total = max_signature(cib) + 1
println("Basin analysis (total = $total starting scenarios)")
println()
println("Fixed points (sorted by basin size):")

perm = sortperm(basins, rev=true)
for (rank, idx) in enumerate(perm)
    u = fixed_points[idx]
    pct = round(100.0 * basins[idx] / total, digits=1)
    println("  #$rank  sig=$(signature(cib, u))  basin=$(basins[idx])  ($pct%)")
    for (i, desc) in enumerate(cib.descriptors)
        println("       $desc = $(cib.variants[desc][u[i] + 1])")
    end
end

if cycle_count > 0
    pct = round(100.0 * cycle_count / total, digits=1)
    println()
    println("Scenarios trapped in non-fixed-point cycles: $cycle_count ($pct%)")
end

"""
End-to-end CIB workflow on the small `CIB_global` example (3 descriptors,
10 variants, 36-scenario space).

Run from the repo root:
    julia --project=. examples/01_global.jl
"""

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files", "CIB_global.scw")

cib = load_scw(SAMPLE)

println("Descriptors and variants")
for desc in cib.descriptors
    println("  $desc: ", join(cib.variants[desc], ", "))
end
println()

println("Total scenario space size: ", max_signature(cib) + 1)
println("Consistent scenarios (fixed points): ", length(cib.consistentScenarios))

for (k, u) in enumerate(cib.consistentScenarios)
    println()
    println("  #$k  (signature $(signature(cib, u)))")
    for (i, desc) in enumerate(cib.descriptors)
        println("    $desc = $(cib.variants[desc][u[i] + 1])")
    end
end

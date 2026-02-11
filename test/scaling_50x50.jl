"""
Test scaling: how does Julia perform with increasing MC threshold
on the 50x50 CIM? Determines whether exhaustive search is practical.

Run with:
    julia --project=. test/scaling_50x50.jl
"""

using CrossImpactBalances

const scw_path = joinpath(@__DIR__, "sample_files", "bench_50x50.scw")

println("="^70)
println("50x50 CIM: MC threshold scaling test")
println("  Scenario space: 60,466,176")
println("="^70)

# Warmup
_ = load_scw(scw_path; mc_threshold=1000)

for mc in [10_000, 50_000, 100_000, 500_000, 1_000_000]
    t0 = time_ns()
    cib = load_scw(scw_path; mc_threshold=mc)
    elapsed = (time_ns() - t0) / 1e9

    sort!(cib.kernel, by=u->signature(cib, u))
    sigs = [signature(cib, u) for u in cib.kernel]

    # Verify all are true fixed points
    all_valid = all(u -> CrossImpactBalances.succession_step(cib, u) == u, cib.kernel)

    println("\n  mc_threshold = $(lpad(string(mc), 10, ' '))")
    println("    Time:    $(round(elapsed, digits=3))s")
    println("    Kernel:  $(length(cib.kernel)) fixed points")
    println("    Sigs:    $sigs")
    println("    Valid:   $all_valid")
    println("    Rate:    $(round(mc / elapsed, digits=0)) scenarios/sec")
end

println("\n", "="^70)
println("Extrapolated exhaustive time: see rate × 60,466,176")
println("="^70)

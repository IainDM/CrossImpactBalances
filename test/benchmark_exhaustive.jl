"""
Exhaustive search benchmark for the 50x50 CIM.

Enumerates all 60,466,176 scenarios to find every fixed point — no sampling,
no uncertainty. Uses all available threads for parallelism.

Run with:
    julia -t auto --project=. test/benchmark_exhaustive.jl
    julia -t 8   --project=. test/benchmark_exhaustive.jl
"""

using CrossImpactBalances

const scw_path = joinpath(@__DIR__, "sample_files", "bench_50x50.scw")

println("="^70)
println("Exhaustive search: 50x50 CIM")
println("  Scenario space: 60,466,176")
println("  Threads:        $(Threads.nthreads())")
println("="^70)

# ── Warmup (compilation) ────────────────────────────────────────────────────
println("\nWarmup (compilation)...")
_ = load_scw(joinpath(@__DIR__, "sample_files", "CIB_global.scw"))

# ── Exhaustive search ───────────────────────────────────────────────────────
println("\nExhaustive search ($(Threads.nthreads()) threads)...")
t0 = time_ns()
cib_ex = load_scw(scw_path)
t_ex = (time_ns() - t0) / 1e9

sort!(cib_ex.consistentScenarios, by=u->signature(cib_ex, u))
ex_sigs = [signature(cib_ex, u) for u in cib_ex.consistentScenarios]

# Verify every result is a true fixed point
all_valid = all(u -> CrossImpactBalances.succession_step(cib_ex, u) == u, cib_ex.consistentScenarios)

println("  Time:    $(round(t_ex, digits=2))s")
println("  Found:   $(length(cib_ex.consistentScenarios)) fixed points (DEFINITIVE)")
println("  Sigs:    $ex_sigs")
println("  Valid:   $all_valid")
println("  Rate:    $(round(60_466_176 / t_ex, digits=0)) scenarios/sec")

# ── Basin analysis ──────────────────────────────────────────────────────────
println("\nBasin analysis (cached chain-following)...")
t0 = time_ns()
fps, basins, cyc = find_basins(load_scw(scw_path; kernel=Vector{Vector{Int}}()))
t_ba = (time_ns() - t0) / 1e9

ba_sigs = [signature(cib_ex, u) for u in fps]
perm = sortperm(basins, rev=true)

println("  Time:    $(round(t_ba, digits=2))s")
println("  Fixed points:   $(length(fps))")
println("  Cycle scenarios: $cyc")
println("  Basins:")
for i in perm
    pct = round(100 * basins[i] / 60_466_176, digits=2)
    println("    sig $(ba_sigs[i]):  $(basins[i]) scenarios ($pct%)")
end
println("  Check: basins + cycles = $(sum(basins) + cyc) (expect 60,466,176)")

# ── Summary ─────────────────────────────────────────────────────────────────
println("\n", "-"^70)
println("SUMMARY")
println("-"^70)
println("  Exhaustive search:  $(length(ex_sigs)) fixed points in $(round(t_ex, digits=2))s")
println("  Basin analysis:     $(length(fps)) fixed points + basins in $(round(t_ba, digits=2))s")
println("  All verified:       $all_valid")
println("\n", "="^70)

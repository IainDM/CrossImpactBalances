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
_ = load_scw(joinpath(@__DIR__, "sample_files", "CIB_global.scw"); exhaustive=true)

# ── MC sampling baseline (10K) ──────────────────────────────────────────────
println("\nMC sampling baseline (mc_threshold=10,000)...")
t0 = time_ns()
cib_mc = load_scw(scw_path; mc_threshold=10_000)
t_mc = (time_ns() - t0) / 1e9

sort!(cib_mc.kernel, by=u->signature(cib_mc, u))
mc_sigs = [signature(cib_mc, u) for u in cib_mc.kernel]
println("  Time:    $(round(t_mc, digits=3))s")
println("  Found:   $(length(cib_mc.kernel)) fixed points")
println("  Sigs:    $mc_sigs")

# ── Exhaustive search ───────────────────────────────────────────────────────
println("\nExhaustive search ($(Threads.nthreads()) threads)...")
t0 = time_ns()
cib_ex = load_scw(scw_path; exhaustive=true)
t_ex = (time_ns() - t0) / 1e9

sort!(cib_ex.kernel, by=u->signature(cib_ex, u))
ex_sigs = [signature(cib_ex, u) for u in cib_ex.kernel]

# Verify every result is a true fixed point
all_valid = all(u -> CrossImpactBalances.succession_step(cib_ex, u) == u, cib_ex.kernel)

println("  Time:    $(round(t_ex, digits=2))s")
println("  Found:   $(length(cib_ex.kernel)) fixed points (DEFINITIVE)")
println("  Sigs:    $ex_sigs")
println("  Valid:   $all_valid")
println("  Rate:    $(round(60_466_176 / t_ex, digits=0)) scenarios/sec")

# ── Summary ─────────────────────────────────────────────────────────────────
missed = setdiff(Set(ex_sigs), Set(mc_sigs))
println("\n", "-"^70)
println("SUMMARY")
println("-"^70)
println("  MC sampling (10K):  $(length(mc_sigs)) fixed points in $(round(t_mc, digits=3))s")
println("  Exhaustive search:  $(length(ex_sigs)) fixed points in $(round(t_ex, digits=2))s")
if !isempty(missed)
    println("  MC missed:          $(length(missed)) fixed points (sigs: $(sort(collect(missed))))")
end
println("  All verified:       $all_valid")
println("\n", "="^70)

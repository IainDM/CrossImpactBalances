"""
Benchmark: 50x50 CIM (15 descriptors, ~60M scenarios).
Compares Julia vs Python CIBSA. Since both use MC sampling with
independent RNGs, we verify fixed-point validity rather than
requiring exact kernel agreement.

Run with:
    julia --project=. test/benchmark_50x50.jl
"""

using CrossImpactBalances
using Test

# ── JSON parser (shared with mcp/server.jl, test/benchmark.jl) ──────────────

include(joinpath(@__DIR__, "..", "mcp", "json.jl"))

# ── Main ─────────────────────────────────────────────────────────────────────

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")
const scw_path = joinpath(SAMPLE_DIR, "bench_50x50.scw")
const expected = parse_json_file(joinpath(SAMPLE_DIR, "bench_50x50_expected.json"))

nvars = Int[v for v in expected["nvariants"]]
total = expected["total_scenarios"]
py_kernel = [Int[v for v in row] for row in expected["kernel"]]
py_sigs = sort(Int[v for v in expected["kernel_sigs"]])
py_time_find = expected["python_time_find_consistent_s"]

println("="^70)
println("50x50 CIM Benchmark: Julia vs Python CIBSA")
println("  $(length(nvars)) descriptors, variants=$nvars")
println("  CIM: $(sum(nvars))×$(sum(nvars)), scenario space: $total")
println("="^70)

@testset "50x50 benchmark" begin

    # ── Warmup (compilation) ──
    println("\nWarmup run (includes compilation)...")
    t0 = time_ns()
    cib_warmup = load_scw(scw_path)
    t_warmup = (time_ns() - t0) / 1e9
    println("  Warmup: $(round(t_warmup, digits=3))s, $(length(cib_warmup.consistentScenarios)) consistent")

    # ── Benchmark: find_consistent (3 runs, take median) ──
    println("\nBenchmark runs (3 trials)...")
    times_find = Float64[]
    local cib
    for trial in 1:3
        t0 = time_ns()
        cib = load_scw(scw_path)
        elapsed = (time_ns() - t0) / 1e9
        push!(times_find, elapsed)
        println("  Trial $trial: $(round(elapsed, digits=4))s, $(length(cib.consistentScenarios)) consistent")
    end
    sort!(times_find)
    jl_time_find = times_find[2]  # median

    # Sort kernel
    sort!(cib.consistentScenarios, by=u->signature(cib, u))
    jl_sigs = [signature(cib, u) for u in cib.consistentScenarios]

    # ── Verify Julia's kernel: every entry is a true fixed point ──
    println("\nVerifying Julia kernel ($(length(cib.consistentScenarios)) scenarios)...")
    @testset "Julia fixed points" begin
        for u in cib.consistentScenarios
            v = CrossImpactBalances.succession_step(cib, u)
            @test v == u
        end
    end

    # ── Verify Python's kernel: every entry is a fixed point under Julia's CIM ──
    println("Verifying Python kernel ($(length(py_kernel)) scenarios)...")
    @testset "Python fixed points verified in Julia" begin
        for u_py in py_kernel
            # Python uses 0-indexed variants, Julia uses 0-indexed too internally
            v = CrossImpactBalances.succession_step(cib, u_py)
            @test v == u_py
        end
    end

    # ── Check kernel overlap ──
    jl_sig_set = Set(jl_sigs)
    py_sig_set = Set(py_sigs)
    overlap = intersect(jl_sig_set, py_sig_set)
    println("\n  Julia kernel:  $(length(jl_sigs)) scenarios, sigs=$jl_sigs")
    println("  Python kernel: $(length(py_sigs)) scenarios, sigs=$py_sigs")
    println("  Overlap:       $(length(overlap)) scenarios")

    # Both should find at least one fixed point
    @test length(cib.consistentScenarios) >= 1
    @test length(py_kernel) >= 1

    # ── Timing report ──
    speedup_find = py_time_find / max(jl_time_find, 1e-9)

    println("\n", "-"^60)
    println("TIMING RESULTS")
    println("-"^60)
    println()
    println("  find_consistent ($(sum(nvars))×$(sum(nvars)) CIM, $(total) scenario space):")
    println("    Python:  $(round(py_time_find, digits=2))s")
    println("    Julia:   $(round(jl_time_find, digits=4))s")
    println("    Speedup: $(round(speedup_find, digits=1))x")
end

println("\n", "="^70)
println("Done.")

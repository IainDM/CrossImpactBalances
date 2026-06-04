"""
Benchmark: load random .scw files, find consistent scenarios, and compare
against Python CIBSA expected results. Measures and reports timing.

Run with:
    julia --project=. test/benchmark.jl
"""

using CrossImpactBalances
using Test

# ── Minimal JSON parser (shared with mcp/server.jl) ─────────────────────────

include(joinpath(@__DIR__, "..", "mcp", "json.jl"))

# ── Main benchmark ──────────────────────────────────────────────────────────

const SAMPLE_DIR = joinpath(@__DIR__, "sample_files")
const EXPECTED = parse_json_file(joinpath(SAMPLE_DIR, "benchmark_expected.json"))

println("="^70)
println("CrossImpactBalances.jl vs Python CIBSA — Benchmark")
println("="^70)

@testset "Benchmarks" begin

for (name, expected) in sort(collect(EXPECTED), by=x->x[2]["total_scenarios"])
    scw_path = joinpath(SAMPLE_DIR, "$name.scw")
    nvars = Int[v for v in expected["nvariants"]]
    total = expected["total_scenarios"]
    n_kernel_expected = expected["n_kernel"]
    expected_sigs = sort(Int[v for v in expected["kernel_sigs"]])
    expected_ipm = expected["inner_product_matrix"]
    py_time_find = expected["python_time_find_consistent_s"]
    py_time_ipm = expected["python_time_ipm_s"]

    println("\n", "-"^60)
    println("$name: $(length(nvars)) descriptors, variants=$nvars, total=$total")

    @testset "$name" begin
        # ── Warmup run (includes compilation) ──
        cib_warmup = load_scw(scw_path)

        # ── Benchmark: find_consistent (3 runs, take median) ──
        times_find = Float64[]
        local cib
        for trial in 1:3
            t0 = time_ns()
            cib = load_scw(scw_path)
            push!(times_find, (time_ns() - t0) / 1e9)
        end
        sort!(times_find)
        jl_time_find = times_find[2]  # median

        # Sort kernel by signature for deterministic ordering
        sort!(cib.kernel, by=u->signature(cib, u))
        kernel_sigs = [signature(cib, u) for u in cib.kernel]

        # ── Verify agreement ──
        @test length(cib.kernel) == n_kernel_expected
        @test kernel_sigs == expected_sigs

        println("  Consistent scenarios: $(length(cib.kernel))")
        println("  Kernel sigs: $kernel_sigs")

        # ── Benchmark: inner_product_matrix ──
        times_ipm = Float64[]
        local M
        for trial in 1:3
            t0 = time_ns()
            M = inner_product_matrix(cib)
            push!(times_ipm, (time_ns() - t0) / 1e9)
        end
        sort!(times_ipm)
        jl_time_ipm = times_ipm[2]

        # ── Verify IPM agreement ──
        if n_kernel_expected > 0
            for (i, row) in enumerate(expected_ipm)
                for (j, val) in enumerate(row)
                    @test M[i, j] == val
                end
            end
        end

        # ── Report timing ──
        speedup_find = py_time_find / max(jl_time_find, 1e-9)
        speedup_ipm = py_time_ipm / max(jl_time_ipm, 1e-9)

        println()
        println("  find_consistent:")
        println("    Python:  $(round(py_time_find, digits=4))s")
        println("    Julia:   $(round(jl_time_find, digits=6))s")
        println("    Speedup: $(round(speedup_find, digits=1))x")
        println()
        println("  inner_product_matrix:")
        println("    Python:  $(round(py_time_ipm, digits=6))s")
        println("    Julia:   $(round(jl_time_ipm, digits=6))s")
        println("    Speedup: $(round(speedup_ipm, digits=1))x")
    end
end

end  # @testset

println("\n", "="^70)
println("Done.")

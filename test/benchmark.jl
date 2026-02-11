"""
Benchmark: load random .scw files, find consistent scenarios, and compare
against Python CIBSA expected results. Measures and reports timing.

Run with:
    julia --project=. test/benchmark.jl
"""

using CrossImpactBalances
using Test

# ── Minimal JSON parser (avoids needing JSON.jl) ────────────────────────────

function parse_json_file(path::String)
    text = read(path, String)
    return parse_json_value(text, 1)[1]
end

function skip_ws(s, i)
    while i <= length(s) && s[i] in (' ', '\t', '\n', '\r')
        i += 1
    end
    return i
end

function parse_json_value(s, i)
    i = skip_ws(s, i)
    i > length(s) && error("Unexpected end of JSON")
    c = s[i]
    if c == '"'
        return parse_json_string(s, i)
    elseif c == '{'
        return parse_json_object(s, i)
    elseif c == '['
        return parse_json_array(s, i)
    elseif c == 't'
        return (true, i + 4)
    elseif c == 'f'
        return (false, i + 5)
    elseif c == 'n'
        return (nothing, i + 4)
    elseif isdigit(c) || c == '-'
        return parse_json_number(s, i)
    else
        error("Unexpected character at position $i: '$(s[i])' ($(Int(s[i])))")
    end
end

function parse_json_string(s, i)
    i += 1  # skip opening "
    buf = IOBuffer()
    while i <= length(s) && s[i] != '"'
        if s[i] == '\\'
            i += 1
            if s[i] == 'n'; write(buf, '\n')
            elseif s[i] == 't'; write(buf, '\t')
            elseif s[i] == '"'; write(buf, '"')
            elseif s[i] == '\\'; write(buf, '\\')
            else write(buf, s[i])
            end
        else
            write(buf, s[i])
        end
        i += 1
    end
    return (String(take!(buf)), i + 1)  # skip closing "
end

function parse_json_number(s, i)
    j = i
    # Optional leading minus
    if j <= length(s) && s[j] == '-'
        j += 1
    end
    # Integer part
    while j <= length(s) && isdigit(s[j])
        j += 1
    end
    is_float = false
    # Fractional part
    if j <= length(s) && s[j] == '.'
        is_float = true
        j += 1
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    # Exponent part
    if j <= length(s) && s[j] in ('e', 'E')
        is_float = true
        j += 1
        if j <= length(s) && s[j] in ('+', '-')
            j += 1
        end
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    numstr = s[i:j-1]
    val = is_float ? parse(Float64, numstr) : parse(Int, numstr)
    return (val, j)
end

function parse_json_array(s, i)
    i += 1  # skip [
    arr = Any[]
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == ']'
        return (arr, i + 1)
    end
    while true
        val, i = parse_json_value(s, i)
        push!(arr, val)
        i = skip_ws(s, i)
        if i > length(s) || s[i] == ']'
            return (arr, i + 1)
        end
        i += 1  # skip comma
    end
end

function parse_json_object(s, i)
    i += 1  # skip {
    obj = Dict{String, Any}()
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == '}'
        return (obj, i + 1)
    end
    while true
        i = skip_ws(s, i)
        key, i = parse_json_string(s, i)
        i = skip_ws(s, i)
        i += 1  # skip colon
        val, i = parse_json_value(s, i)
        obj[key] = val
        i = skip_ws(s, i)
        if i > length(s) || s[i] == '}'
            return (obj, i + 1)
        end
        i += 1  # skip comma
    end
end

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

"""
Standalone analysis script invoked by the Node.js MCP wrapper.

Usage:
    julia --startup-file=no --project=<project_dir> mcp/run_analysis.jl <file_path>

Prints the formatted analysis result to stdout.
Exits with code 1 and prints to stderr on failure.
"""

using CrossImpactBalances

function format_results(cib, fps, basins, cycle_count, total)
    buf = IOBuffer()

    println(buf, "Cross-Impact Balance Analysis")
    println(buf, "=" ^ 50)
    println(buf)

    # Problem size
    println(buf, "Descriptors: $(cib.ndesc)")
    print(buf, "Variants per descriptor: [")
    print(buf, join(cib.nvariants, ", "))
    println(buf, "]")
    println(buf, "Total scenario space: $total")
    println(buf)

    # Descriptor listing
    println(buf, "Descriptors and variants:")
    for desc in cib.descriptors
        vars = cib.variants[desc]
        println(buf, "  $desc: $(join(vars, ", "))")
    end
    println(buf)

    # Fixed points sorted by basin size (descending)
    perm = sortperm(basins, rev=true)
    println(buf, "Consistent Scenarios (Fixed Points): $(length(fps))")
    println(buf, "-" ^ 50)

    for (rank, idx) in enumerate(perm)
        fp = fps[idx]
        basin = basins[idx]
        pct = round(100.0 * basin / total; digits=2)

        println(buf)
        println(buf, "Scenario #$rank (basin: $basin scenarios, $pct%)")

        for (j, desc) in enumerate(cib.descriptors)
            variant_name = cib.variants[desc][fp[j] + 1]  # 0-based to 1-based
            println(buf, "  $desc = $variant_name")
        end
    end

    # Cycles
    if cycle_count > 0
        println(buf)
        cyc_pct = round(100.0 * cycle_count / total; digits=2)
        println(buf, "Scenarios in cycles (non-convergent): $cycle_count ($cyc_pct%)")
    end

    # Summary
    println(buf)
    converging = sum(basins)
    println(buf, "Coverage: $converging converging + $cycle_count cycling = $total total")

    return String(take!(buf))
end

function main()
    if length(ARGS) < 1
        println(stderr, "Usage: run_analysis.jl <file_path>")
        exit(1)
    end

    file_path = ARGS[1]

    if !isfile(file_path)
        println(stderr, "File not found: $file_path")
        exit(1)
    end

    cib = load_scw(file_path; kernel=Vector{Vector{Int}}())
    fps, basins, cycle_count = find_basins(cib)
    total = max_signature(cib) + 1

    text = format_results(cib, fps, basins, cycle_count, total)
    print(text)
end

main()

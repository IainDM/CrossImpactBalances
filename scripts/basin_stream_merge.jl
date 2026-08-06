# Combine basin_stream_worker.jl partials into the full exact result.
#
#   julia --project=. scripts/basin_stream_merge.jl model.scw part1.tsv part2.tsv ...
#
# Validates that the partials' ranges tile the model's whole scenario space
# 0:(scenario_count-1) exactly — no gap, no overlap, no partial from a
# different model — then sums the per-fixed-point counts and the cycle
# counts, checks the coverage invariant (sum of basins + cycles == scenario
# count), and prints the final table, largest basin first. This is exactly
# the result find_basins would have produced in one piece.

using CrossImpactBalances

function read_partial(path)
    model = ""
    range = (-1, -1)
    counts = Dict{Int,Int}()
    cycles = -1
    for line in eachline(path)
        fields = split(rstrip(line), '\t')
        isempty(fields[1]) && continue
        if fields[1] == "M"
            model = String(fields[2])
        elseif fields[1] == "R"
            range = (parse(Int, fields[2]), parse(Int, fields[3]))
        elseif fields[1] == "F"
            sig = parse(Int, fields[2])
            haskey(counts, sig) && error("$path: duplicate fixed point $sig")
            counts[sig] = parse(Int, fields[3])
        elseif fields[1] == "C"
            cycles = parse(Int, fields[2])
        else
            error("$path: unrecognised line kind $(repr(fields[1]))")
        end
    end
    model == "" && error("$path: missing M line")
    range[1] < 0 && error("$path: missing R line")
    cycles < 0 && error("$path: missing C line")
    return model, range, counts, cycles
end

function main(args)
    if length(args) < 2
        println(stderr, "usage: basin_stream_merge.jl model.scw part1.tsv [part2.tsv ...]")
        return 2
    end
    scwPath = args[1]
    cib = load_scw(scwPath; kernel=Vector{Vector{Int}}())
    totalScenarios = scenario_count(cib)
    lastValid = max_signature(cib)

    partials = [read_partial(path) for path in args[2:end]]

    # Every partial must be from this model, and the ranges must tile
    # 0:lastValid exactly.
    expectedModel = basename(scwPath)
    for (index, (model, range, _, _)) in enumerate(partials)
        model == expectedModel || error(
            "$(args[1 + index]) is a partial of $(repr(model)), not $(repr(expectedModel))")
        (0 <= range[1] <= range[2] <= lastValid) || error(
            "$(args[1 + index]): range $(range) outside 0:$(lastValid)")
    end
    ranges = sort!([partial[2] for partial in partials])
    ranges[1][1] == 0 || error("the space is not covered: no partial starts at 0")
    for index in 2:length(ranges)
        expectedNext = ranges[index - 1][2] + 1
        ranges[index][1] == expectedNext || error(
            ranges[index][1] > expectedNext ?
            "gap: no partial covers $(expectedNext):$(ranges[index][1] - 1)" :
            "overlap: $(ranges[index]) begins inside $(ranges[index - 1])")
    end
    ranges[end][2] == lastValid || error(
        "the space is not covered: nothing beyond $(ranges[end][2]) (need $(lastValid))")

    merged = Dict{Int,Int}()
    totalCycles = 0
    for (_, _, counts, cycles) in partials
        for (sig, count) in counts
            merged[sig] = get(merged, sig, 0) + count
        end
        totalCycles += cycles
    end

    covered = sum(values(merged); init=0) + totalCycles
    Int128(covered) == totalScenarios || error(
        "coverage invariant failed: $(covered) != $(totalScenarios)")

    println("# exact basins of ", expectedModel, " — merged from ", length(partials),
            " partials")
    println("# total scenarios: ", totalScenarios)
    println("# consistent scenarios: ", length(merged))
    println("# starts ending in cycles: ", totalCycles)
    println("rank\tsignature\tbasin_size\tbasin_pct\tscenario")
    bySize = sort!(collect(merged); by=pair -> -pair[2])
    for (rank, (sig, size)) in enumerate(bySize)
        scenario = inv_signature(cib, sig)
        pct = round(100 * size / Float64(totalScenarios); sigdigits=4)
        println(rank, '\t', sig, '\t', size, '\t', pct, '\t', scenario)
    end
    return 0
end

exit(main(ARGS))

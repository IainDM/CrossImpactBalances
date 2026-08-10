#!/usr/bin/env julia
# Verify JuCIB against Wolfgang Weimer-Jehle's ScenarioWizard corpus.
#
# Three blocks, in this order:
#   1. INVENTORY  — every .scw: shape, scenario count, whether a solution set
#                   ships with it.
#   2. EXPORTS    — every .sl: does its header's declared count match the rows
#                   it actually holds, and is every listed scenario a genuine
#                   fixed point of the matrix it is paired with?
#   3. RE-SEARCH  — JuCIB finds the kernel itself and compares, reporting the
#                   two directions SEPARATELY: `missing` (ScenarioWizard found
#                   it, we did not — that would mean find_consistent is
#                   unsound, the most serious result this script can produce)
#                   and `extra` (we found it, ScenarioWizard did not — either
#                   our bug or a genuine finding about ScenarioWizard). A single
#                   "verified" boolean hides the difference; do not add one.
#
# The corpus is NOT in the repository — see test/WWJ_CORPUS.md for what it is,
# why, and how to stage a copy at test/wwj_corpus/.
#
# Run:  julia -t auto --project=. test/verify_wwj.jl [name ...]
#       (no arguments = every model; names = just those, e.g. `D60 B80`)
#
# Exit codes:  0 all good · 1 a discrepancy · 2 no corpus staged

using CrossImpactBalances
const CIBmod = CrossImpactBalances

const CORPUS = joinpath(@__DIR__, "wwj_corpus")
const OUTPUT = joinpath(@__DIR__, "bench_results_wwj.json")

# .sl exports that do NOT belong to the .scw of the same name: right descriptor
# count and variant shape, but not fixed points of that matrix — evidently
# exported from an earlier revision of it. Their disagreement is the expected
# result, so it must not fail the run. See test/WWJ_CORPUS.md § Known-bad files.
const STALE_SL = ["N45", "N50"]

# Some kernels are too large to hold as a list — B90's is 1.9e7 scenarios and
# B100's is 1.9e10 — but their influence maps split, so the size can be
# computed exactly as a product over islands and never built. Take that route
# whenever it is available AND nothing needs the scenarios themselves:
#
#   paired with a .sl  →  enumerate, because the comparison needs the actual
#                         set (the largest is B80's 18,432 — no trouble)
#   one island only    →  enumerate, because "decomposing" would be the same
#                         search with extra steps
#   otherwise          →  product over islands
#
# A memory decision, not a capability one: the answer is exact either way.
decompose(cib, hasReference) = !hasReference && length(influence_structure(cib).components) > 1

if !isdir(CORPUS)
    println("""
        No corpus staged at $(CORPUS).

        These are Wolfgang Weimer-Jehle's own ScenarioWizard files and they are
        not redistributed with this repository. See test/WWJ_CORPUS.md for what
        the corpus contains and how to put a copy in place.""")
    exit(2)
end

# ─── Minimal JSON writing (same shape as test/verify_sw.jl) ──────────────────
jstr(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
jval(x::Bool)           = x ? "true" : "false"
jval(x::Integer)        = string(x)
jval(x::AbstractFloat)  = (isnan(x) || isinf(x)) ? "null" : string(x)
jval(::Nothing)         = "null"
jval(x::AbstractString) = jstr(x)
jval(x::Vector{<:Pair}) = "{" * join([jstr(k) * ":" * jval(v) for (k, v) in x], ",") * "}"
jval(x::AbstractVector) = "[" * join([jval(v) for v in x], ",") * "]"

# The count ScenarioWizard wrote into line 3 of the .sl. Comparing it with the
# rows actually present is what catches a truncated export.
function declared_solution_count(slPath)
    for (lineNumber, line) in enumerate(eachline(slPath))
        lineNumber == 3 && return tryparse(Int, strip(line))
    end
    return nothing
end

models = sort([f[1:end-4] for f in readdir(CORPUS) if endswith(f, ".scw")])
if !isempty(ARGS)
    wanted = Set(ARGS)
    unknown = setdiff(wanted, Set(models))
    isempty(unknown) || error("no such model in the corpus: $(join(sort(collect(unknown)), ", "))")
    models = filter(in(wanted), models)
end

records = Any[]
discrepancies = String[]

println("JuCIB vs ScenarioWizard — Weimer-Jehle corpus")
println("Julia $(VERSION), $(Threads.nthreads()) thread(s), $(length(models)) model(s)\n")

# ── 1. Inventory ─────────────────────────────────────────────────────────────
println("── Inventory " * "─"^66)
println(rpad("model", 8), rpad("desc", 6), rpad("var", 6), rpad("scenarios", 26),
        rpad(">Int64", 8), "solution set")
loaded = Dict{String,CIB}()
for name in models
    cib = load_scw(joinpath(CORPUS, "$name.scw"); kernel = Vector{Vector{Int}}())
    loaded[name] = cib
    total = scenario_count(cib)
    slPath = joinpath(CORPUS, "$name.sl")
    status = !isfile(slPath) ? "none" : (name in STALE_SL ? "STALE (see WWJ_CORPUS.md)" : "paired")
    println(rpad(name, 8), rpad(cib.numberOfDescriptors, 6), rpad(cib.numberOfDimensions, 6),
            rpad(string(total), 26), rpad(total > Int128(typemax(Int)), 8), status)
end

# ── 2. ScenarioWizard's exports ──────────────────────────────────────────────
println("\n── ScenarioWizard's exports " * "─"^51)
println(rpad("model", 8), rpad("declared", 10), rpad("rows", 8), rpad("in range", 10),
        rpad("unique", 8), "fixed points")
solutions = Dict{String,Vector{Vector{Int}}}()
for name in models
    slPath = joinpath(CORPUS, "$name.sl")
    isfile(slPath) || continue
    cib = loaded[name]
    sols = load_solutions(cib, slPath)
    solutions[name] = sols

    declared = declared_solution_count(slPath)
    inRange = all(u -> length(u) == cib.numberOfDescriptors &&
                       all(0 .<= u .< cib.numberOfVariants), sols)
    unique = allunique(sols)
    # Only meaningful once the indices are in range — an out-of-range scenario
    # cannot be scored at all.
    fixedPoints = inRange ? count(u -> CIBmod.succession_step(cib, u) == u, sols) : 0

    note = if !inRange
        "n/a (indices out of range)"
    elseif fixedPoints == length(sols)
        "$(fixedPoints)/$(length(sols)) ✓"
    else
        "$(fixedPoints)/$(length(sols)) " * (name in STALE_SL ? "(stale, expected)" : "*** NOT ALL ***")
    end
    if fixedPoints != length(sols) && !(name in STALE_SL)
        push!(discrepancies, "$name: only $fixedPoints of $(length(sols)) exported scenarios are fixed points")
    end
    if declared != length(sols)
        push!(discrepancies, "$name.sl: header declares $declared solutions but holds $(length(sols))")
    end
    println(rpad(name, 8), rpad(something(declared, "?"), 10), rpad(length(sols), 8),
            rpad(inRange, 10), rpad(unique, 8), note)
end

# ── 3. JuCIB's own search ────────────────────────────────────────────────────
println("\n── JuCIB's own search " * "─"^57)
println(rpad("model", 8), rpad("kernel", 12), rpad("SW", 8), rpad("missing", 9),
        rpad("extra", 8), rpad("time", 9), "how")
for name in models
    cib = loaded[name]
    reference = get(solutions, name, nothing)
    isStale = name in STALE_SL

    useIslands = decompose(cib, reference !== nothing)
    local kernelSize, kernel, how
    elapsed = @elapsed begin
        if useIslands
            # Exact, but as a product: the islands are independent, so the
            # kernel is their Cartesian product and we never build it.
            islands = split_cib(cib)
            islandKernels = [length(find_consistent(island)) for island in islands]
            kernelSize = prod(Int128.(islandKernels))
            kernel = nothing
            how = "split_cib ($(length(islands)) islands, core $(maximum(islandKernels)))"
        else
            kernel = find_consistent(cib)
            kernelSize = Int128(length(kernel))
            how = "exhaustive"
        end
    end

    missingCount, extraCount = nothing, nothing
    if reference !== nothing && kernel !== nothing
        found = Set(kernel)
        expected = Set(reference)
        missingCount = count(u -> !(u in found), reference)
        extraCount = count(u -> !(u in expected), kernel)
        if !isStale && (missingCount > 0 || extraCount > 0)
            push!(discrepancies,
                  "$name: $missingCount ScenarioWizard scenarios not found, $extraCount extra found")
            for u in Iterators.take(Iterators.filter(v -> !(v in found), reference), 3)
                push!(discrepancies, "    missing: $u")
            end
            for u in Iterators.take(Iterators.filter(v -> !(v in expected), kernel), 3)
                push!(discrepancies, "    extra:   $u")
            end
        end
        isStale && (how *= " (vs stale export — mismatch expected)")
    end

    push!(records, Pair{String,Any}[
        "file" => name,
        "n_descriptors" => cib.numberOfDescriptors,
        "n_variants" => cib.numberOfDimensions,
        "scenario_count" => string(scenario_count(cib)),
        "over_int64" => scenario_count(cib) > Int128(typemax(Int)),
        "sw_declared" => reference === nothing ? nothing :
                         declared_solution_count(joinpath(CORPUS, "$name.sl")),
        "sw_rows" => reference === nothing ? nothing : length(reference),
        "sw_stale" => isStale,
        "kernel_size" => string(kernelSize),
        "missing_count" => missingCount,
        "extra_count" => extraCount,
        "method" => useIslands ? "split_cib" : "exhaustive",
        "time_s" => round(elapsed; digits = 2),
        "threads" => Threads.nthreads(),
        "julia_version" => string(VERSION),
    ])

    println(rpad(name, 8), rpad(string(kernelSize), 12),
            rpad(reference === nothing ? "—" : length(reference), 8),
            rpad(something(missingCount, "—"), 9), rpad(something(extraCount, "—"), 8),
            rpad(string(round(elapsed; digits = 2)) * "s", 9), how)
    flush(stdout)
end

open(OUTPUT, "w") do io
    write(io, jval(records))
end
println("\nWrote $(relpath(OUTPUT, dirname(@__DIR__)))  ($(length(records)) models)")

if isempty(discrepancies)
    println("\nAll good: no missing scenarios, no extra scenarios, no truncated exports.")
    exit(0)
else
    println("\n*** $(length(discrepancies)) discrepancy line(s):")
    for line in discrepancies
        println("  ", line)
    end
    exit(1)
end

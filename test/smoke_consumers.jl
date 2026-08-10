# ══ Smoke tests for code OUTSIDE src/ that depends on the engine's API ═════
#
# Everything else in this suite tests `src/`. But several things ship with the
# package and call into it from outside: the desktop app, and the documented
# examples. Nothing else here touches them, so when the engine's API changes
# they can sit broken indefinitely — and the first sign of it is a failed
# release build, or a user.
#
# That is not hypothetical. Two real cases motivated this file:
#
#   * `app/src/CIBApp.jl` kept using `cib.ndesc` after that field was renamed
#     to `numberOfDescriptors`. Precompiling the app does not execute the line,
#     so CI stayed green; the MSI would have installed and then thrown on the
#     user's first click.
#
#   * The quick-start example in README.md and docs/src/index.md iterated
#     `cib.kernel`, also renamed. That is the first code a new user copies.
#
#   * `capi/src/CIBCApi.jl` reported the scenario space as `max_signature + 1`.
#     That THROWS past typemax(Int64), and the throw came back from `cib_load`,
#     so a model that size could not be loaded through the C API at all — and
#     `cib_consistent`, which handles those models perfectly well, was
#     unreachable for them. Nothing here executed a line of that module: it is
#     compiled into libcib by PackageCompiler and driven over ctypes from
#     `python/crossimpactbalances/_native.py`.
#
# So these tests deliberately EXECUTE things rather than just loading them.
# They assert only on shape and internal consistency, not exact numbers —
# the point is "does it still run and produce something coherent", because
# the engine's actual answers are pinned thoroughly elsewhere.

using Test
using CrossImpactBalances

const REPO = normpath(joinpath(@__DIR__, ".."))

# Run `f()` with stdout silenced — the examples and app are chatty and their
# output would drown the test report. Errors still propagate normally.
function quietly(f)
    original = stdout
    (rd, wr) = redirect_stdout()
    reader = @async read(rd, String)   # drain the pipe so a chatty script can't block
    try
        return f()
    finally
        redirect_stdout(original)
        close(wr)
        wait(reader)
        close(rd)
    end
end

"""Write a FORCING CHAIN of `n` binary descriptors as a ScenarioWizard `.scw`.

The same matrix as `forcing_chain` in property_tests.jl — descriptor `i` gives
+1 to descriptor `i+1`'s same variant — so the kernel is
`{all-zeros, all-ones}` whatever `n` is, and the pruning is total. Written to a
file because the C API loads from a path, and a model of 2^70 scenarios is not
something to ship as a fixture.
"""
function write_forcing_chain_scw(path::String, n::Int)
    cim = zeros(Int, 2n, 2n)
    for i in 1:n-1, v in 0:1
        cim[2(i - 1) + v + 1, 2i + v + 1] = 1
    end
    open(path, "w") do io
        println(io, "\$ ScenarioWizard 4.0")
        println(io, basename(path))
        for i in 1:n
            println(io, "&D$i")
            println(io, " -V$(i)a")
            println(io, " -V$(i)b")
        end
        # One '#' closes the header; four more step past the metadata sections
        # the parser skips, leaving it reading the matrix.
        for _ in 1:5
            println(io, "#")
        end
        for row in 1:2n
            println(io, join(cim[row, :], ","))
        end
    end
    return path
end

@testset "Desktop app (app/src/CIBApp.jl)" begin
    # CIBApp is a separate package, but its only non-stdlib dependency is
    # CrossImpactBalances itself, so it loads straight from source here
    # without needing app/'s own environment instantiated.
    sandbox = Module(:CIBAppSmoke)
    Base.include(sandbox, joinpath(REPO, "app", "src", "CIBApp.jl"))
    App = getfield(sandbox, :CIBApp)

    samples = sort(filter(f -> endswith(f, ".scw"),
                          readdir(joinpath(REPO, "test", "sample_files"))))
    @test !isempty(samples)

    # Only the small models — this is a smoke test, not a benchmark. The
    # 400-million-scenario files are covered by the engine's own tests.
    small = filter(samples) do f
        cib = load_scw(joinpath(REPO, "test", "sample_files", f);
                       kernel = Vector{Vector{Int}}())
        # scenario_count, not max_signature + 1: max_signature THROWS past
        # typemax(Int64), so a big model dropped into sample_files would break
        # this filter rather than be filtered out by it. Identical value for
        # every file here today.
        scenario_count(cib) <= 100_000
    end
    @test !isempty(small)

    for f in small
        path = joinpath(REPO, "test", "sample_files", f)
        @testset "$f" begin
            cib = App.load_from_text(read(path, String))
            state = App.AppState()

            # ── Find Consistent Scenarios ──
            consistent = quietly(() -> App.analyze_consistent(cib))
            @test consistent["mode"] == "consistent"
            @test consistent["total"] == max_signature(cib) + 1
            @test consistent["count"] == length(consistent["scenarios"])
            @test consistent["descriptors"] == cib.descriptors
            # The app must agree with the engine it wraps.
            @test consistent["count"] == length(find_consistent(cib))
            # variant_names must render one NAME per descriptor — this is the
            # call that silently broke on the cib.ndesc rename.
            for s in consistent["scenarios"]
                @test length(s["variants"]) == cib.numberOfDescriptors
                @test all(v -> v isa AbstractString && !isempty(v), s["variants"])
            end
            @test App.to_json(consistent) isa AbstractString

            # ── Find Basins ──
            basins = quietly(() -> App.analyze_basins(cib, state))
            @test basins["mode"] == "basins"
            @test basins["count"] == consistent["count"]        # same attractors
            @test basins["covered"] + basins["cycles"] == basins["total"]
            @test issorted([s["basin_size"] for s in basins["scenarios"]]; rev = true)
            for s in basins["scenarios"]
                @test length(s["variants"]) == cib.numberOfDescriptors
            end
            @test App.to_json(basins) isa AbstractString

            # ── CSV export ── (populated as a side effect of analyze_basins)
            @test !isempty(state.last_csv)
            @test state.last_csv_name == "basin_analysis.csv"
            rows = split(strip(state.last_csv), '\n')
            header = findfirst(r -> startswith(r, "rank,"), rows)
            @test header !== nothing
            # One data row per consistent scenario, each with the right arity.
            @test length(rows) - header == basins["count"]
            expected_columns = 2 + cib.numberOfDescriptors + 2
            for r in rows[header:end]
                @test count(==(','), r) == expected_columns - 1
            end

            # ── Estimate Basin Shares ──
            estimate = quietly(() -> App.analyze_estimate(cib, state; samples = 5_000))
            @test estimate["mode"] == "estimate"
            @test estimate["samples"] == 5_000
            # The kernel is exact, so the estimate reports the same attractors
            # the exact analysis found — zero-hit ones included.
            @test estimate["count"] == basins["count"]
            @test sum(s["hits"] for s in estimate["scenarios"]; init = 0) +
                  estimate["cycle_hits"] == estimate["samples"]
            for s in estimate["scenarios"]
                @test length(s["variants"]) == cib.numberOfDescriptors
                @test 0.0 <= s["ci_lo_pct"] <= s["ci_hi_pct"] <= 100.0
            end
            @test App.to_json(estimate) isa AbstractString
            @test state.last_csv_name == "basin_share_estimate.csv"
            @test occursin("share_pct", state.last_csv)

            # ── Transitions (the lever map) ──
            transitions = quietly(() -> App.analyze_transitions(cib, state))
            @test transitions["mode"] == "transitions"
            @test transitions["count"] == basins["count"]      # same consistent scenarios
            @test length(transitions["variants"]) == cib.numberOfDescriptors
            @test transitions["world_index"] == 0              # none supplied
            @test occursin("digraph", transitions["dot"])
            for node in transitions["nodes"]
                @test length(node["variants"]) == cib.numberOfDescriptors
                @test node["kind"] in ("attractor", "cycle", "world")
            end
            for edge in transitions["edges"]
                @test 1 <= edge["from"] <= length(transitions["nodes"])
                @test 1 <= edge["to"] <= length(transitions["nodes"])
                @test !edge["baseline"]                        # no world ⇒ no baseline edge
                for change in edge["changes"]
                    @test change["descriptor"] in cib.descriptors
                    @test change["from"] != change["to"]
                end
            end
            @test App.to_json(transitions) isa AbstractString
            @test state.last_csv_name == "transition_graph.csv"

            # With a world state: one extra node, and a baseline edge from it.
            world = zeros(Int, cib.numberOfDescriptors)
            withWorld = quietly(() -> App.analyze_transitions(cib, state; world = world))
            @test withWorld["world_index"] == length(withWorld["nodes"])
            worldNode = withWorld["nodes"][withWorld["world_index"]]
            @test worldNode["kind"] in ("world", "attractor")   # unless it IS consistent
            if worldNode["kind"] == "world"
                @test any(e -> e["baseline"] && e["from"] == withWorld["world_index"],
                          withWorld["edges"])
            end

            # ── Structure ──
            structure = quietly(() -> App.analyze_structure(cib))
            @test structure["mode"] == "structure"
            @test !isempty(structure["islands"])
            for island in structure["islands"]
                @test !isempty(island["descriptors"])
            end
            # Descriptors partition across islands, none lost or repeated.
            all_members = reduce(vcat, [island["descriptors"]
                                        for island in structure["islands"]])
            @test sort(all_members) == sort(cib.descriptors)
            @test App.to_json(structure) isa AbstractString
        end
    end

    # A model far past the table method's memory must come back as guidance
    # (the basins_too_big shape the page renders), never as a crash — and the
    # estimator must work on the very same model. Exercised through the same
    # dispatch the /analyze route uses.
    @testset "huge model: guidance + estimate" begin
        nv = vcat(fill(4, 11), fill(3, 13))          # 2^22 · 3^13 ≈ 6.7e12 scenarios
        ndim = sum(nv)
        offsets = cumsum(vcat(0, nv[1:end-1]))
        cim = zeros(Int, ndim, ndim)
        cim[1, offsets[2] + 1] = 1                    # one real influence, else silent
        # A pre-supplied kernel (the all-zeros scenario is a fixed point of
        # this nearly-silent matrix): estimate_basins would otherwise run
        # find_consistent, and this degenerate matrix has ~10^12 consistent
        # scenarios — a hazard of the synthetic test model, not of real ones.
        huge = CIB(["D$i" for i in 1:24],
                   Dict("D$i" => ["V$(i)_$j" for j in 1:nv[i]] for i in 1:24),
                   nv, cim, permutedims(cim), ndim, 24,
                   [zeros(Int, 24)], offsets)
        state = App.AppState()
        @test_throws ArgumentError App.analyze_basins(huge, state)

        est = quietly(() -> App.analyze_estimate(huge, state; samples = 2_000))
        @test est["mode"] == "estimate"
        @test est["total"] == 6_687_075_336_192      # < 2^53: crosses as a number
        @test sum(s["hits"] for s in est["scenarios"]; init = 0) +
              est["cycle_hits"] == 2_000

        st = quietly(() -> App.analyze_structure(huge))
        @test st["mode"] == "structure"
        @test length(st["islands"]) > 1              # the silent descriptors split off
        @test App.to_json(st) isa AbstractString

        # The lever map costs the same at any scenario-space size, so it needs
        # no guidance path at all — it simply runs. (The kernel here comes from
        # the model's own consistentScenarios, set above; find_consistent must
        # never be let loose on this degenerate matrix.)
        tr = quietly(() -> App.analyze_transitions(huge, state; world = ones(Int, 24)))
        @test tr["mode"] == "transitions"
        @test tr["total"] == 6_687_075_336_192
        @test tr["world_index"] == length(tr["nodes"])
        @test !isempty(tr["edges"])
        @test App.to_json(tr) isa AbstractString
    end

    # A kernel too large to draw must come back as guidance from every mode
    # that lists one consistent scenario per row — carrying the EXACT size,
    # which resolve_kernel reads off the influence map rather than by building
    # the kernel. Without this the page tries to render every row: the WWJ
    # corpus has models with 1.2 million consistent scenarios (~470 MB of JSON)
    # and one with 19,327,352,832, which no machine can hold at all.
    @testset "kernel too big: guidance, not a 4096-row table" begin
        # A silent matrix: nothing influences anything, so every one of the
        # 2^12 scenarios is consistent AND the model splits into 12 islands —
        # the decomposing route, where the size is a product over islands.
        n = 12
        nv = fill(2, n)
        cim = zeros(Int, 2n, 2n)
        wide = CIB(["D$i" for i in 1:n],
                   Dict("D$i" => ["V$(i)_$j" for j in 1:2] for i in 1:n),
                   nv, cim, permutedims(cim), 2n, n,
                   Vector{Vector{Int}}(), cumsum(vcat(0, nv[1:end-1])))
        state = App.AppState()
        @test length(influence_structure(wide).components) == n
        @test length(find_consistent(wide)) == 4096       # what is NOT rendered

        for (mode, call) in (("consistent", () -> App.analyze_consistent(wide)),
                             ("estimate",
                              () -> App.analyze_estimate(wide, state; samples = 2_000)),
                             ("transitions", () -> App.analyze_transitions(wide, state)))
            refusal = quietly(call)
            @test refusal["mode"] == "kernel_too_big"
            @test refusal["requested"] == mode
            @test refusal["kernel"] == 4096
            @test refusal["total"] == 4096
            @test refusal["cap"] == App.KERNEL_DISPLAY_CAP
            @test !isempty(refusal["message"])
            # The refusal is what the page renders, so it must encode — and
            # stay small, which is the whole point of not sending the kernel.
            json = App.to_json(refusal)
            @test json isa AbstractString
            @test length(json) < 10_000
        end

        # Exact Basins is the one mode that reaches the cap holding a finished
        # result, so between the two caps it gives up the table but keeps the
        # export — the whole analysis, in the one place that can hold it.
        drawn = App.AppState()
        withCsv = quietly(() -> App.analyze_basins(wide, drawn))
        @test withCsv["mode"] == "kernel_too_big"
        @test withCsv["requested"] == "basins"
        @test withCsv["kernel"] == 4096
        @test withCsv["export"] === true
        @test occursin("Export CSV", withCsv["message"])
        @test occursin("\"export\":true", App.to_json(withCsv))
        @test drawn.last_csv_name == "basin_analysis.csv"
        @test count(==('\n'), drawn.last_csv) == 4096 + 5   # 4 comments, header, rows

        # Past the export cap not even that — just the size, and nothing stale
        # left behind for /export.csv to hand over in its place.
        bare = App.AppState()
        noCsv = quietly(() -> App.analyze_basins(wide, bare; export_cap = 100))
        @test noCsv["mode"] == "kernel_too_big"
        @test noCsv["kernel"] == 4096
        @test noCsv["export"] === false
        @test isempty(bare.last_csv)

        # Structure is the route the message points at, so it must still run.
        @test quietly(() -> App.analyze_structure(wide))["mode"] == "structure"

        # A model carrying its own kernel is taken at its word, never searched:
        # here that kernel is one scenario, so the guard lets it straight
        # through even though a search would have found 4096.
        carried = CIB(wide.descriptors, wide.variants, wide.numberOfVariants,
                      wide.cim, wide.cim_t, 2n, n,
                      [zeros(Int, n)], wide.desc_offsets)
        @test quietly(() -> App.analyze_consistent(carried))["count"] == 1
    end
end

@testset "C API (capi/src/CIBCApi.jl)" begin
    # CIBCApi is a separate package whose only non-stdlib dependency is
    # CrossImpactBalances, so it loads straight from source — no compiled
    # library, no PackageCompiler, no Python.
    sandbox = Module(:CApiSmoke)
    Base.include(sandbox, joinpath(REPO, "capi", "src", "CIBCApi.jl"))
    CApi = getfield(sandbox, :CIBCApi)

    # Call a cib_* function the way ctypes does: NUL-terminated C strings in,
    # a malloc'd C string out (freed here), decoded with the module's own
    # parser — so what these tests see is exactly what crosses the boundary.
    function call_c(f, args...)
        held = String[]
        cargs = Any[]
        for arg in args
            if arg isa AbstractString
                push!(held, String(arg) * "\0")
                push!(cargs, Cstring(pointer(held[end])))
            else
                push!(cargs, arg)
            end
        end
        reply = GC.@preserve held f(cargs...)
        text = unsafe_string(reply)
        CApi.cib_free_string(reply)
        return CApi.parse_json_value(text, 1)[1]
    end

    as_json(u) = "[" * join(u, ",") * "]"
    as_indices(x) = Int[Int(v) for v in x]

    @testset "small model: the C surface agrees with the engine it wraps" begin
        path = joinpath(REPO, "test", "sample_files", "CIB_global.scw")
        cib = load_scw(path; kernel = Vector{Vector{Int}}())

        loaded = call_c(CApi.cib_load, path)
        @test loaded["ok"]
        handle = Cint(loaded["handle"])
        @test loaded["descriptors"] == cib.descriptors
        @test loaded["n_descriptors"] == cib.numberOfDescriptors
        # Comfortably under 2^53, so the count crosses as a JSON number.
        @test loaded["n_scenarios"] isa Integer
        @test loaded["n_scenarios"] == scenario_count(cib) == 36

        kernel = call_c(CApi.cib_consistent, handle, "{}")
        @test kernel["ok"]
        found = [as_indices(u) for u in kernel["scenarios"]]
        @test found == find_consistent(cib)
        @test [signature(cib, u) for u in found] == [13, 16, 20, 21]

        for u in found
            sig = call_c(CApi.cib_signature, handle, as_json(u))
            @test sig["signature"] == signature(cib, u)
            back = call_c(CApi.cib_inv_signature, handle, string(signature(cib, u)))
            @test as_indices(back["scenario"]) == u
        end

        # From a consistent scenario: one step, and it is its own successor.
        atRest = call_c(CApi.cib_succession, handle,
                        "{\"scenario\":" * as_json(found[1]) * "}")
        @test atRest["ok"]
        @test atRest["cycle_length"] == 1
        @test [as_indices(s) for s in atRest["steps"]] == [found[1]]

        # From all-first-variants this model enters a TWO-cycle. Reporting a
        # cycle as length 1 is what `converged` is built from downstream, so
        # the length is the assertion that matters here.
        walk = call_c(CApi.cib_succession, handle, "{\"scenario\":[0,0,0]}")
        @test walk["ok"]
        @test [as_indices(s) for s in walk["steps"]] ==
              [[0, 0, 0], [1, 1, 2], [0, 1, 1], [1, 1, 2]]
        @test walk["cycle_length"] == 2

        basins = call_c(CApi.cib_basins, handle, "{}")
        @test basins["ok"]
        @test basins["total"] == 36

        CApi.cib_free(handle)
    end

    # THE REGRESSION. 70 binary descriptors is 1.18e21 scenarios: past
    # typemax(Int64), so there are no signatures, and every line of this
    # surface that reached for one used to be either a refusal or a silent
    # wrap. find_consistent searches such a model anyway — branch-and-bound
    # numbers no scenario — and the C API has to be able to ask it.
    @testset "past typemax(Int64)" begin
        n = 70
        path = write_forcing_chain_scw(joinpath(mktempdir(), "forcing_chain.scw"), n)
        cib = load_scw(path; kernel = Vector{Vector{Int}}())
        expectedCount = Int128(2)^n
        @test scenario_count(cib) == expectedCount
        @test expectedCount > Int128(typemax(Int))
        @test_throws ArgumentError max_signature(cib)     # what cib_load used to call

        loaded = call_c(CApi.cib_load, path)
        @test loaded["ok"]                                # it used to come back false
        handle = Cint(loaded["handle"])
        @test loaded["n_descriptors"] == n
        # Past 2^53 the count crosses as a decimal string, exact to the digit:
        # a consumer that parses JSON numbers into doubles would otherwise
        # round the last four digits of this one away.
        @test loaded["n_scenarios"] isa AbstractString
        @test parse(Int128, loaded["n_scenarios"]) == expectedCount

        kernel = call_c(CApi.cib_consistent, handle, "{}")
        @test kernel["ok"]
        @test [as_indices(u) for u in kernel["scenarios"]] == [zeros(Int, n), ones(Int, n)]

        # Signatures are the Int128 ones. For the all-ones scenario the true
        # answer is 2^70 - 1; the Int64 `signature` wraps it to exactly -1, and
        # that is the number this surface used to hand back.
        @test signature(cib, ones(Int, n)) == -1
        highest = call_c(CApi.cib_signature, handle, as_json(ones(Int, n)))
        @test highest["signature"] isa AbstractString
        @test parse(Int128, highest["signature"]) == expectedCount - 1
        lowest = call_c(CApi.cib_signature, handle, as_json(zeros(Int, n)))
        @test lowest["signature"] == 0

        # Inverting is Int64 work whatever the model size, so it refuses by
        # name past that — and still answers for any signature that fits.
        refusal = call_c(CApi.cib_inv_signature, handle, "\"$(expectedCount - 1)\"")
        @test refusal["ok"] == false
        @test occursin("past typemax(Int64)", refusal["error"])
        inRange = call_c(CApi.cib_inv_signature, handle, "12345")
        @test inRange["ok"]
        @test as_indices(inRange["scenario"]) == inv_signature(cib, 12345)

        # Succession must reach the attractor and report cycle_length 1. Keyed
        # on Int64 signatures it could report a cycle that is not there, two
        # different scenarios having wrapped onto one key.
        start = [isodd(i) ? 1 : 0 for i in 1:n]
        walk = call_c(CApi.cib_succession, handle,
                      "{\"scenario\":" * as_json(start) * "}")
        @test walk["ok"]
        @test walk["cycle_length"] == 1
        @test as_indices(walk["steps"][end]) in (zeros(Int, n), ones(Int, n))

        # find_basins genuinely cannot run at this size. Its refusal is what
        # must come back — not a wrapped total from the reply builder.
        basins = call_c(CApi.cib_basins, handle, "{}")
        @test basins["ok"] == false
        @test occursin("estimate_basins", basins["error"])

        CApi.cib_free(handle)
    end
end

@testset "Documented examples run" begin
    examples = sort(filter(f -> endswith(f, ".jl"), readdir(joinpath(REPO, "examples"))))
    @test !isempty(examples)
    for f in examples
        @testset "$f" begin
            # A fresh module each time, so examples cannot collide over the
            # `const SAMPLE` and struct definitions they each declare.
            sandbox = Module(Symbol("Example_", replace(f, r"\W" => "_")))
            quietly(() -> Base.include(sandbox, joinpath(REPO, "examples", f)))
            @test true      # reaching here means it ran without throwing
        end
    end
end

@testset "Quick-start snippet in README.md and docs/src/index.md" begin
    # Extract the `cib.<field>` accesses and function calls the documented
    # snippets use, and check each still exists. This catches the exact drift
    # that left `cib.kernel` in both files after the field was renamed.
    cib = load_scw(joinpath(REPO, "test", "sample_files", "CIB_global.scw"))

    for doc in (joinpath(REPO, "README.md"), joinpath(REPO, "docs", "src", "index.md"))
        @testset "$(basename(doc))" begin
            text = read(doc, String)
            # Julia code fences only. `\r?` because these files are a mix of
            # CRLF and LF line endings, and matching only \n silently finds
            # nothing rather than failing usefully.
            blocks = [m.captures[1] for m in eachmatch(r"```julia\r?\n(.*?)```"s, text)]
            @test !isempty(blocks)

            fields = Set{Symbol}()
            for b in blocks, m in eachmatch(r"\bcib\.([A-Za-z_][A-Za-z0-9_]*)", b)
                push!(fields, Symbol(m.captures[1]))
            end
            @test !isempty(fields)
            for name in fields
                @test hasfield(CIB, name) || hasproperty(cib, name)
            end

            # Any exported name the snippets call must still be exported.
            exported = Set(names(CrossImpactBalances))
            for b in blocks, m in eachmatch(r"\b([a-z_][A-Za-z0-9_!]*)\(", b)
                fn = Symbol(m.captures[1])
                fn in exported && @test isdefined(CrossImpactBalances, fn)
            end
        end
    end
end

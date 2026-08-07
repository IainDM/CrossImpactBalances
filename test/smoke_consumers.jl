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
        max_signature(cib) + 1 <= 100_000
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

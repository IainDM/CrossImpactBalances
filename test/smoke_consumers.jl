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
            rows = split(strip(state.last_csv), '\n')
            header = findfirst(r -> startswith(r, "rank,"), rows)
            @test header !== nothing
            # One data row per consistent scenario, each with the right arity.
            @test length(rows) - header == basins["count"]
            expected_columns = 2 + cib.numberOfDescriptors + 2
            for r in rows[header:end]
                @test count(==(','), r) == expected_columns - 1
            end
        end
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

# Wolfgang Weimer-Jehle's ScenarioWizard corpus: 38 real models from the author
# of the method, and the solution sets ScenarioWizard itself produced for 19 of
# them. Provenance, the full table and the three known-bad exports:
# test/WWJ_CORPUS.md
#
# THE CORPUS IS NOT IN THE REPOSITORY. It is Weimer-Jehle's material and is not
# redistributed, so test/wwj_corpus/ is gitignored and everything below skips
# when it is absent. That is why this file guards on `isdir` rather than
# asserting the data exists.
#
# ── NEVER USE `signature` IN THIS FILE ──────────────────────────────────────
# Fourteen of these models have more scenarios than typemax(Int64), including
# seven of the nineteen paired ones. `signature` accumulates in Int and wraps
# silently: on B80 its 18,432 distinct scenarios collapse onto 9,024 distinct
# values. A test written the way the sample-file round-trip at
# test/runtests.jl:250 is written — comparing sorted signature vectors — would
# compare wrapped garbage against wrapped garbage and pass.
#
# Compare the scenarios themselves. `sort` on Vector{Vector{Int}} is
# lexicographic and total, so `sort(a) == sort(b)` is exact set equality at any
# model size, and `Set` comparison works too.
#
# ── WHY SOME OF IT IS BEHIND AN ENV VAR ─────────────────────────────────────
# `Pkg.test()` runs Julia with `--check-bounds=yes`, which disables every
# `@inbounds` in the search kernels. That costs the exhaustive searches here
# roughly 7x: N40 takes 4.8 s in an ordinary session and 34.8 s under the test
# runner. Searching all 38 models would add 18 minutes to a single-threaded CI
# leg, so the two slowest paired models, the models with no reference solution,
# and the four whose kernels have to be composed rather than enumerated are
# held back behind JUCIB_WWJ_FULL=1.
#
# What stays on by default is the claim itself: 17 of the 19 paired models
# re-searched and compared set-for-set, plus every cheap check on all 38. About
# 1m55s here. `test/verify_wwj.jl` does the complete job, all 38 models and both
# disagreement directions, outside the suite and without the bounds penalty.
const WWJ_FULL = get(ENV, "JUCIB_WWJ_FULL", "") in ("1", "true", "yes")

const WWJ_DIR = joinpath(@__DIR__, "wwj_corpus")

# The corpus, ordered by size as test/WWJ_CORPUS.md lists it.
#   kernel = consistent scenarios JuCIB finds (pinned; for the four whose
#            kernels are too large to hold, the exact product over islands)
#   nsol   = solutions in the paired .sl, or `nothing` for the 19 models with
#            no usable reference — including N45 and N50, whose .sl files are
#            stale exports and are checked separately below
const WWJ_MODELS = [
    (name="N20",   ndesc=20,  nvar=54,  kernel=33,             nsol=nothing),
    (name="N25",   ndesc=25,  nvar=67,  kernel=3,              nsol=3),
    (name="N25a",  ndesc=25,  nvar=67,  kernel=3,              nsol=nothing),
    (name="N25b",  ndesc=25,  nvar=67,  kernel=2,              nsol=nothing),
    (name="N25c",  ndesc=25,  nvar=67,  kernel=10,             nsol=nothing),
    (name="D30a",  ndesc=30,  nvar=75,  kernel=13,             nsol=13),
    (name="N30",   ndesc=30,  nvar=81,  kernel=15,             nsol=15),
    (name="D35",   ndesc=35,  nvar=88,  kernel=51,             nsol=51),
    (name="D35a",  ndesc=35,  nvar=88,  kernel=15,             nsol=nothing),
    (name="N35",   ndesc=35,  nvar=94,  kernel=3,              nsol=3),
    (name="D40",   ndesc=40,  nvar=100, kernel=110,            nsol=110),
    (name="D40_a", ndesc=40,  nvar=120, kernel=379,            nsol=nothing),
    (name="D40a",  ndesc=40,  nvar=100, kernel=3240,           nsol=nothing),
    (name="N40",   ndesc=40,  nvar=108, kernel=15,             nsol=15),
    (name="N40a",  ndesc=40,  nvar=108, kernel=18,             nsol=18),
    (name="D45",   ndesc=45,  nvar=113, kernel=86,             nsol=86),
    (name="D45a",  ndesc=45,  nvar=113, kernel=12,             nsol=nothing),
    (name="N45",   ndesc=45,  nvar=121, kernel=62,             nsol=nothing),
    (name="N45a",  ndesc=45,  nvar=121, kernel=2,              nsol=nothing),
    (name="N45b",  ndesc=45,  nvar=121, kernel=18,             nsol=nothing),
    (name="N45c",  ndesc=45,  nvar=121, kernel=42,             nsol=nothing),
    (name="B50",   ndesc=50,  nvar=100, kernel=172,            nsol=172),
    (name="D50",   ndesc=50,  nvar=125, kernel=470,            nsol=470),
    (name="D50a",  ndesc=50,  nvar=125, kernel=35,             nsol=nothing),
    (name="N50",   ndesc=50,  nvar=135, kernel=11,             nsol=nothing),
    (name="D55",   ndesc=55,  nvar=138, kernel=6340,           nsol=nothing),
    (name="D55a",  ndesc=55,  nvar=138, kernel=140,            nsol=140),
    (name="D55b",  ndesc=55,  nvar=138, kernel=328,            nsol=328),
    (name="D55c",  ndesc=55,  nvar=138, kernel=620,            nsol=620),
    (name="B60",   ndesc=60,  nvar=120, kernel=928,            nsol=928),
    (name="D60",   ndesc=60,  nvar=150, kernel=194,            nsol=194),
    (name="D60b",  ndesc=60,  nvar=120, kernel=382,            nsol=382),
    (name="B70",   ndesc=70,  nvar=140, kernel=3156,           nsol=3156),
    (name="B80",   ndesc=80,  nvar=160, kernel=18432,          nsol=18432),
    (name="B90",   ndesc=90,  nvar=180, kernel=18874368,       nsol=nothing),
    (name="B100",  ndesc=100, nvar=200, kernel=19327352832,    nsol=nothing),
    (name="B100b", ndesc=100, nvar=200, kernel=1195648,        nsol=nothing),
    (name="B100c", ndesc=100, nvar=200, kernel=26992,          nsol=nothing),
]

# Kernels too large to hold as a list — B90's is 1.9e7 scenarios and B100's is
# 1.9e10 — but their influence maps split, so the size is exactly the product
# over islands and is computed without ever building it. B100b joins them for
# the same reason (its core island alone has 597,824). B100c has one island, so
# there is nothing to decompose and it is searched whole.
const WWJ_BY_ISLANDS = ["B90", "B100", "B100b"]

# The two paired models whose re-search does not fit a default CI leg. Measured
# under `--check-bounds=yes` at one thread: N40a 72.8 s, D55b ~100 s, against
# 30 s or less for every other pair. Both run under JUCIB_WWJ_FULL=1 and in
# test/verify_wwj.jl; ScenarioWizard's answers for them are still checked to be
# genuine fixed points by default, which is the soundness half.
const WWJ_SLOW_PAIRS = ["N40a", "D55b"]

# .sl exports that do not belong to the .scw of the same name — see
# test/WWJ_CORPUS.md § Known-bad files. Pinned two-sided: `fixed` is how many of
# their scenarios really are fixed points of the matrix they are named after,
# and it is 0 for both. A change that made these start passing would be as
# suspicious as one that made a good pair start failing.
const WWJ_STALE_SL = [(name="N45", rows=29, fixed=0), (name="N50", rows=93, fixed=0)]

# Models whose scenario space passes typemax(Int64), so they have no signatures:
# find_consistent still searches them, max_signature and find_basins must not.
const WWJ_OVER_INT64 = ["B100", "B100b", "B100c", "B70", "B80", "B90",
                        "D50", "D50a", "D55", "D55a", "D55b", "D55c", "D60", "N50"]

wwj_path(name, ext) = joinpath(WWJ_DIR, name * ext)
wwj_load(name) = load_scw(wwj_path(name, ".scw"); kernel = Vector{Vector{Int}}())

# The count ScenarioWizard wrote into line 3 of the .sl. Comparing it against
# the rows actually present is what catches a truncated export — it is exactly
# the check that identifies B90.sl as unusable.
function wwj_declared_count(slPath)
    for (lineNumber, line) in enumerate(eachline(slPath))
        lineNumber == 3 && return tryparse(Int, strip(line))
    end
    return nothing
end

if !isdir(WWJ_DIR)
    @info """Weimer-Jehle corpus not staged — skipping its tests.
             These are third-party files that this repository does not
             redistribute; see test/WWJ_CORPUS.md to put a copy in place."""
else
@testset "Weimer-Jehle ScenarioWizard corpus" begin

    @testset "manifest integrity" begin
        # Every row has a file and every file has a row: catches a lost fixture
        # and an unregistered addition with the same assertion.
        @test length(WWJ_MODELS) == 38
        onDisk = Set(f[1:end-4] for f in readdir(WWJ_DIR) if endswith(f, ".scw"))
        @test Set(m.name for m in WWJ_MODELS) == onDisk

        expectedSl = Set(m.name for m in WWJ_MODELS if m.nsol !== nothing)
        union!(expectedSl, Set(s.name for s in WWJ_STALE_SL))
        @test Set(f[1:end-3] for f in readdir(WWJ_DIR) if endswith(f, ".sl")) == expectedSl

        # B90.sl is deliberately excluded: 22 descriptors wide, so it belongs to
        # no model here, and truncated at 100,000 of a declared 1,241,136.
        @test !isfile(wwj_path("B90", ".sl"))
    end

    @testset "every .scw parses with the right shape" begin
        for model in WWJ_MODELS
            cib = wwj_load(model.name)
            @test cib.numberOfDescriptors == model.ndesc
            @test length(cib.numberOfVariants) == model.ndesc
            @test cib.numberOfDimensions == model.nvar
            @test sum(cib.numberOfVariants) == model.nvar
            @test size(cib.cim) == (model.nvar, model.nvar)
            @test cib.cim_t == permutedims(cib.cim)
            @test cib.desc_offsets == cumsum(vcat(0, cib.numberOfVariants[1:end-1]))
            @test scenario_count(cib) == prod(Int128.(cib.numberOfVariants))
            # Verifies the WWJ_OVER_INT64 pin rather than trusting it.
            @test (scenario_count(cib) > Int128(typemax(Int))) == (model.name in WWJ_OVER_INT64)
        end
    end

    @testset "models past typemax(Int64) refuse the signature-driven analyses" begin
        # The real-matrix twin of the synthetic guard in large_space_tests.jl.
        for name in WWJ_OVER_INT64
            cib = wwj_load(name)
            @test_throws ArgumentError max_signature(cib)
            @test_throws ArgumentError find_basins(cib)
            @test_throws ArgumentError find_consistent(cib; algorithm=:sweep)
            @test scenario_count(cib) > Int128(typemax(Int))   # but counting is exact
        end
    end

    @testset "ScenarioWizard's solutions are genuine fixed points" begin
        # The cheap direction, and it covers all 25,136 reference scenarios —
        # 18,432 of them B80's — in well under a second.
        for model in WWJ_MODELS
            model.nsol === nothing && continue
            slPath = wwj_path(model.name, ".sl")
            cib = load_scw(wwj_path(model.name, ".scw"); sl_file = slPath)
            sols = cib.consistentScenarios

            @test length(sols) == model.nsol
            # The header's own count must match what the file holds.
            @test wwj_declared_count(slPath) == length(sols)
            @test allunique(sols)
            @test all(u -> length(u) == cib.numberOfDescriptors, sols)
            @test all(u -> all(0 .<= u .< cib.numberOfVariants), sols)

            # One assertion per model rather than per scenario: 25,136 @test
            # records would swamp the report, and findfirst keeps a failure
            # diagnosable.
            offender = findfirst(u -> CrossImpactBalances.succession_step(cib, u) != u, sols)
            @test offender === nothing
        end
    end

    @testset "JuCIB finds exactly ScenarioWizard's kernel" begin
        # The expensive direction, and the point of the whole corpus: nothing
        # ScenarioWizard found that JuCIB misses, and nothing JuCIB finds that
        # ScenarioWizard missed. Spaces up to 1.2e24.
        #
        # algorithm=:bnb, not :auto. On a budget trip :auto falls back to the
        # odometer sweep, which on B50 would mean walking 1.1e15 scenarios —
        # a CI hang rather than a test failure. :bnb has no fallback.
        for model in WWJ_MODELS
            model.nsol === nothing && continue
            (!WWJ_FULL && model.name in WWJ_SLOW_PAIRS) && continue
            cib = wwj_load(model.name)
            kernel = find_consistent(cib; algorithm = :bnb)
            reference = load_solutions(cib, wwj_path(model.name, ".sl"))
            @test length(kernel) == model.kernel
            @test sort(kernel) == sort(reference)      # exact set equality
        end
        WWJ_FULL || @info "WWJ: re-searched 17 of 19 paired models; " *
                          "$(join(WWJ_SLOW_PAIRS, ", ")) need JUCIB_WWJ_FULL=1 " *
                          "(or test/verify_wwj.jl)."
    end

    @testset "stale ScenarioWizard exports (N45.sl, N50.sl)" begin
        for stale in WWJ_STALE_SL
            cib = wwj_load(stale.name)
            sols = load_solutions(cib, wwj_path(stale.name, ".sl"))

            # They load cleanly — right width, indices in range — which is
            # exactly why the mismatch is not obvious without checking.
            @test length(sols) == stale.rows
            @test all(u -> length(u) == cib.numberOfDescriptors, sols)
            @test all(u -> all(0 .<= u .< cib.numberOfVariants), sols)

            @test count(u -> CrossImpactBalances.succession_step(cib, u) == u, sols) == stale.fixed
        end
    end

    # ── From here on: JUCIB_WWJ_FULL=1 only ─────────────────────────────────
    # These prove nothing about ScenarioWizard parity — there is no reference
    # solution to compare against — so they are the right things to drop from a
    # default CI leg. test/verify_wwj.jl always runs them.
    if !WWJ_FULL
        @info "WWJ: set JUCIB_WWJ_FULL=1 for the models with no reference " *
              "solution and the composed 10^27/10^30 kernels."
    else

    @testset "kernels for the models with no reference solution" begin
        # No ScenarioWizard answer to compare against, so these are JuCIB's own
        # result, pinned as regression values. Every scenario is still checked
        # to be a genuine fixed point — the one claim that stands alone.
        for model in WWJ_MODELS
            model.nsol === nothing || continue
            model.name in WWJ_BY_ISLANDS && continue
            cib = wwj_load(model.name)
            kernel = find_consistent(cib; algorithm = :bnb)
            @test length(kernel) == model.kernel
            offender = findfirst(u -> CrossImpactBalances.succession_step(cib, u) != u, kernel)
            @test offender === nothing
        end
    end

    @testset "the 10^27 and 10^30 models, through split_cib" begin
        # B90's kernel is 18.9 million scenarios and B100's is 19.3 billion:
        # exact, but not something to hold in memory (a direct search for B90's
        # passed 15 GB before being killed). Their influence maps split, so the
        # size is the product over islands and never gets built. Do NOT replace
        # this with find_consistent — that is the whole point.
        for model in WWJ_MODELS
            model.name in WWJ_BY_ISLANDS || continue
            cib = wwj_load(model.name)
            islands = split_cib(cib)
            @test length(islands) > 1
            product = prod(Int128(length(find_consistent(island))) for island in islands)
            @test product == Int128(model.kernel)
        end

        # B100c is the control: same 10^30 space, one island, so it is searched
        # whole and lands on a kernel small enough to hold.
        b100c = wwj_load("B100c")
        @test length(influence_structure(b100c).components) == 1
        @test length(find_consistent(b100c; algorithm = :bnb)) == 26992
    end

    end   # WWJ_FULL

    @testset "blank variant names" begin
        # The corpus's name generator runs off the end of the alphabet and then
        # emits nothing at all, so real ScenarioWizard files with empty variant
        # names exist. Index-based addressing is unaffected; name-based
        # addressing silently resolves to the FIRST match instead of raising.
        # Pinned so it cannot change without someone noticing.
        cib = wwj_load("B50")
        @test cib.variants["e"] == ["", ""]
        @test get_impact(cib, "e", "", "f", "") == get_impact(cib, 4, 0, 5, 0)
    end
end
end

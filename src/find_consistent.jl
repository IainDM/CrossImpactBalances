# ══ Finding the consistent scenarios ═══════════════════════════════════════
#
# This file is the entry point and the traffic controller. It decides HOW to
# search; the searching itself lives in sweep.jl and branch_and_bound.jl.
#
# There are two ways to find every consistent scenario, and which one is
# available depends on how the succession rule defines consistency:
#
#   THE SLOW, ALWAYS-CORRECT WAY. Walk every scenario, call the rule's own
#   succession_step, and keep the scenarios that step to themselves. This asks
#   the rule directly, so it works no matter how exotic the rule is. It is
#   also single-threaded and allocates on every scenario.
#
#   THE FAST WAY. Never call succession_step at all. Instead test each
#   scenario directly: "does any variant beat the one in use?" This is a far
#   cheaper question, and it is what lets the sweep and branch-and-bound work.
#   But it is only VALID for rules whose notion of consistency really does
#   reduce to that test — which is what `fixed_point_margin` declares.

"""
    fixed_point_margin(rule) -> Union{Nothing,Int}

Declares that a rule's consistent scenarios can be identified by a simple local test, and how demanding that test is.

Ordinarily a scenario is consistent when no descriptor could do better by switching variant. Declaring a margin `m` relaxes this to "no descriptor could do better *by more than m*". That reflects the fact that real social systems have inertia: they do not reorganise themselves over a trivial improvement. A margin of `0` is the ordinary meaning, and both built-in rules declare it.

The declaration is a promise, not a calculation — the rule is asserting that its consistent scenarios are exactly the scenarios passing that test. Making the promise is what buys the rule the fast searches in [`find_consistent`](@ref). A rule that returns `nothing` (the default for any new rule) makes no promise and is searched the slow, always-correct way.

The name uses "fixed point" — the dynamical-systems term for a consistent scenario, i.e. one whose succession step maps it to itself. The two are the same thing; see the vocabulary note at the top of `CrossImpactBalances.jl`.
"""
fixed_point_margin(::SuccessionRule) = nothing
fixed_point_margin(::GlobalSuccession) = 0
# Sequential succession reaches consistent scenarios by a different route, but
# ends up at exactly the same set as the global rule (provable by induction),
# so it can honestly make the margin-0 promise too.
fixed_point_margin(::SequentialSuccession) = 0

"""
    find_consistent(cib; rule=GlobalSuccession(), algorithm=:auto,
                    bnb_node_budget=nothing) -> Vector{Vector{Int}}

Find every consistent scenario in the model — the model's *kernel*.

`rule` selects the succession dynamics ([`SuccessionRule`](@ref); default [`GlobalSuccession`](@ref)).

`algorithm` selects the search strategy. It only applies to rules that declare a [`fixed_point_margin`](@ref); any other rule is searched one scenario at a time and passing an algorithm is an error.

- `:sweep` — check every scenario, using the incremental odometer walk in `sweep.jl`.
- `:bnb`   — branch-and-bound (`branch_and_bound.jl`): build scenarios one descriptor at a time and discard whole families that provably contain nothing consistent. Exact — it returns the identical answer — and usually looks at a small fraction of the scenarios.
- `:auto`  (default) — the sweep for small models (under 100,000 scenarios, where branch-and-bound's overheads are not repaid), otherwise branch-and-bound with a budget; if the pruning turns out too weak to pay off, the budget trips and the sweep runs instead.

Branch-and-bound never numbers a scenario, so it — and therefore `:auto` above 100,000 scenarios — keeps working past `typemax(Int64)` ≈ 9.2×10¹⁸ scenarios, where signatures stop existing ([`max_signature`](@ref)). Two strategies do not cross that line, and both say so plainly: `:sweep`, which walks signatures in order, and the generic scan used by rules that declare no [`fixed_point_margin`](@ref). Nor is there anything to fall back TO if branch-and-bound exhausts an explicitly-set `bnb_node_budget` there — so `:auto` does not set one: its usual budget of a sixteenth of the scenario count is unreachable at that scale, and the search simply runs to completion.

The result is ordered by ascending signature — decided by comparing the scenarios' variant choices from the last descriptor down, which is the same order and needs no signature to exist. So it does not depend on how many threads happened to be running, at any model size.
"""
function find_consistent(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                         algorithm::Symbol=:auto,
                         bnb_node_budget::Union{Nothing,Int}=nothing)

    return _find_kernel(rule, cib; algorithm=algorithm,
                              bnb_node_budget=bnb_node_budget)
end

"""
    _find_kernel(rule, cib; algorithm, bnb_node_budget)

Choose a search strategy for `rule` and run it. See the notes at the top of this file for the two kinds of search and why not every rule can have the fast one.
"""
function _find_kernel(rule::SuccessionRule, cib::CIB; algorithm::Symbol=:auto,
                            bnb_node_budget::Union{Nothing,Int}=nothing)
    margin = fixed_point_margin(rule)

    # ── No margin declared: the slow, always-correct search ──
    # Walk the scenarios in order, ask the rule where each one goes, and keep
    # the ones that go nowhere. Single-threaded, and it allocates inside
    # succession_step on every scenario, but it makes no assumptions about the
    # rule whatsoever. Neither built-in rule arrives here — both declare 0.
    if margin === nothing
        algorithm === :auto || throw(ArgumentError(
            "algorithm=$(repr(algorithm)) needs a rule with a fixed_point_margin; " *
            "this rule uses the generic scan (algorithm=:auto)"))
        # The scan is signature-driven — one scenario per signature — so it
        # stops where signatures do. Branch-and-bound would carry on past here,
        # but only a rule that declares a margin is allowed to use it.
        totalScenarios = scenario_count(cib)
        totalScenarios <= Int128(typemax(Int)) || throw(ArgumentError(
            "find_consistent: this rule declares no fixed_point_margin, so every " *
            "scenario has to be visited one at a time, and this model has " *
            "$(totalScenarios) scenarios — past typemax(Int64) = $(typemax(Int)) " *
            "there are no signatures to walk. Declaring a fixed_point_margin for " *
            "the rule (see `fixed_point_margin`) unlocks branch-and-bound, which " *
            "needs none; fix_descriptor / influence_structure cut the model down " *
            "instead."))
        consistentScenarios = Vector{Vector{Int}}()
        for currentSignature in 0:max_signature(cib)
            scenario = inv_signature(cib, currentSignature)
            # A scenario that is its own successor is consistent.
            if succession_step(rule, cib, scenario) == scenario
                push!(consistentScenarios, scenario)
            end
        end
        return consistentScenarios
    end

    # ── A margin was declared, so the fast searches are available ──
    algorithm in (:auto, :bnb, :sweep) ||
        throw(ArgumentError("algorithm must be :auto, :bnb or :sweep, got $(repr(algorithm))"))
    # Int128, because branch-and-bound does not number a scenario and so does
    # not stop where signatures do. This count is the only quantity below that
    # has to be representable, and scenario_count is exact to ~1.7e38 (and
    # guards even that). max_signature would refuse outright — and worse, for a
    # model of exactly 2^63 scenarios it returns typemax(Int) WITHOUT refusing,
    # whose `+ 1` wrapped negative, passed the "< 100_000" test below, and sent
    # the search into a sweep of an empty range: an empty kernel, no error.
    numberOfScenarios = scenario_count(cib)

    # Use the sweep when asked for it, and when the model is small enough that
    # branch-and-bound's setup would cost more than it saves.
    if algorithm == :sweep || (algorithm == :auto && numberOfScenarios < 100_000)
        # The sweep is an odometer over signatures, so it needs them to exist.
        # Only an explicit :sweep can fail this test — the :auto arm above is
        # already limited to models of under 100,000 scenarios.
        numberOfScenarios <= Int128(typemax(Int)) || throw(ArgumentError(
            "find_consistent: algorithm=:sweep walks every scenario in signature " *
            "order, and this model has $(numberOfScenarios) scenarios — more than " *
            "typemax(Int64) = $(typemax(Int)), so the signatures it would walk do " *
            "not exist. Use algorithm=:bnb — branch-and-bound needs no signatures, " *
            "and it is what :auto picks at this size anyway — or cut the model down " *
            "with fix_descriptor / influence_structure."))
        # Note the sweep needs only the margin — it never touches the rule.
        return _find_kernel_checkall_fast(cib; margin=margin)
    end

    # Otherwise branch-and-bound, which avoids looking at most scenarios at
    # all by ruling out whole families of them at once.
    #
    # `something(a, b)` returns the first argument that isn't `nothing` — like
    # C#'s ?? null-coalescing operator. An explicit `:bnb` means the caller
    # wants branch-and-bound whatever it costs, so it gets an unlimited
    # budget; `:auto` caps the effort at a sixteenth of the scenario count.
    #
    # The sixteenth is computed in Int128 and then clamped, because past ~1.5e20
    # scenarios it no longer fits an Int. Clamping decides nothing: no search
    # visits typemax(Int) ≈ 9.2e18 nodes, so a budget at or above that was
    # already unreachable in exact arithmetic. The practical consequence is that
    # on a model past Int64 `:auto` runs to completion rather than giving up —
    # which is right, since there is no sweep to fall back to there.
    autoBudget = Int(min(numberOfScenarios ÷ 16, Int128(typemax(Int))))
    nodeBudget = something(bnb_node_budget,
                           algorithm == :bnb ? typemax(Int) : autoBudget)

    # _bnb_search picks the descriptor order, searches in it, and puts the
    # answer back into this model's own descriptor numbering.
    kernel, _ = _bnb_search(cib; node_budget=nodeBudget, margin=margin)

    # Careful: a `nothing` result does NOT mean "this model has no consistent
    # scenarios" — a model with none returns an empty list, perfectly normally.
    # `nothing` means branch-and-bound gave up, having spent its whole node
    # budget without finishing, because on this particular matrix the pruning
    # was too weak to be worth it. (See "why there is a node budget" in
    # branch_and_bound.jl.)
    if !isnothing(kernel)
        return kernel
    end

    # It gave up, so fall back to checking every scenario after all — which
    # only exists while signatures do. Getting here on a bigger model takes an
    # explicitly-set budget: :auto's own budget is unreachable at that scale
    # (see the clamp above), so branch-and-bound runs to completion instead.
    numberOfScenarios <= Int128(typemax(Int)) || throw(ArgumentError(
        "find_consistent: branch-and-bound spent its node budget of $(nodeBudget) " *
        "without finishing, and the sweep it falls back to cannot run on a model of " *
        "$(numberOfScenarios) scenarios (more than typemax(Int64) = $(typemax(Int)), " *
        "so it has no signatures to walk). Raise bnb_node_budget, or leave it unset — " *
        "both :auto and :bnb then let branch-and-bound run to completion — or cut the " *
        "model down with fix_descriptor / influence_structure."))
    return _find_kernel_checkall_fast(cib; margin=margin)
end

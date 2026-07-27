# ══ Succession: how one scenario leads to the next ═════════════════════════
#
# CIB rests on two separate ideas, and it is worth keeping them apart:
#
#   CONSISTENCY  Is this scenario stable? That is a question about a single
#                scenario on its own — does any descriptor wish it had picked
#                a different variant?
#
#   SUCCESSION   If a scenario is NOT stable, where does it go next? That is a
#                question about movement.
#
# The two are independent: "is there somewhere better to be" is a different
# question from "how do you get there". This file is about the second one.
#
# The raw material for both is the impact balance (see impact_balance.jl),
# which scores every variant of every descriptor against the current scenario.
# A descriptor "prefers" whichever of ITS OWN variants scores highest. A
# succession rule decides how those preferences turn into an actual move — and
# the built-in rules differ in exactly one respect: WHEN the scores are
# recalculated.

"""
    SuccessionRule

Succession is one of the two core concepts of CIB — "how do we move from one scenario to the next". The other is consistency — "is there a better scenario than this one". These are separate: whether a better scenario exists is a different question from how you would get there.

To allow flexibility in how succession is calculated, we define an abstract supertype `SuccessionRule`. An abstract type cannot be instantiated itself — it only serves as a parent for concrete rule types, the way an interface or abstract base class would in C# or Python.

To add a new rule, define

    struct MyRule <: SuccessionRule end

(`<:` means "is a subtype of") and a single method

    succession_step(rule::MyRule, cib::CIB, u::Vector{Int}) -> Vector{Int}

and every analysis routine — [`find_consistent`](@ref) and [`find_basins`](@ref) — works with it immediately, because they dispatch on the rule's type. The only assumption is that a rule depends solely on the current scenario, not on the path taken to reach it (in technical terms, that the process is Markov).
"""
abstract type SuccessionRule end

"""
The standard succession rule, as used by ScenarioWizard and the original CIBSA.

Every variant is scored against the *current* scenario, and then every descriptor moves at once to its best-scoring variant. All the moves are decided from the same snapshot of scores, so within a step no descriptor can see what any other descriptor did. Ties favour the current variant, then the lowest variant index.
"""
struct GlobalSuccession <: SuccessionRule end

"""
Sequential (successive / Gauss–Seidel) succession: descriptors move one at a time, in descriptor order, and the scores are recomputed after each move. So descriptor 2 chooses in full knowledge of what descriptor 1 just did, within the same step.

Its *trajectories* differ from [`GlobalSuccession`](@ref) — where an unstable scenario goes next depends on the update order — but its *consistent scenarios* are exactly the same ones. See [`fixed_point_margin`](@ref) for why that matters.
"""
struct SequentialSuccession <: SuccessionRule end

"""
    succession_step(cib, u) -> Vector{Int}
    succession_step(rule, cib, u) -> Vector{Int}

Take one step of succession from scenario `u` and return the scenario it leads to. Omitting `rule` means [`GlobalSuccession`](@ref).

A scenario that steps to itself is a consistent scenario: nobody wants to move, so nothing changes.
"""
function succession_step(::GlobalSuccession, cib::CIB, scenario::Vector{Int})
    # Score every variant ONCE, against the scenario as it stands right now.
    # Every descriptor below then chooses from this single frozen snapshot,
    # which is what makes the moves simultaneous — descriptor 5 cannot tell
    # that descriptor 1 moved, because the scores never change during the loop.
    impactBalance = impact_balance(cib, scenario)
    successor = copy(scenario)

    # impactBalance holds every variant of every descriptor end to end in one
    # long vector, so each descriptor owns a contiguous block of it. We walk
    # the blocks in turn: firstColumn and lastColumn bracket the block of the
    # descriptor currently under consideration, and firstColumn steps past it
    # at the end of each iteration.
    firstColumn = 1
    for descriptorIndex in 1:cib.numberOfDescriptors
        variantCount = cib.numberOfVariants[descriptorIndex]
        lastColumn = firstColumn + variantCount - 1

        # `@view` gives us that block without copying it — like a NumPy view
        # or a C# Span, rather than the fresh array a plain slice allocates.
        descriptorScores = @view impactBalance[firstColumn:lastColumn]

        # Open the contest with the variant the scenario currently uses, so it
        # stands as the incumbent and a challenger has to actually beat it.
        # (+1 because variants are numbered from 0 but arrays index from 1.)
        bestScore = descriptorScores[scenario[descriptorIndex] + 1]

        for variantIndex in 0:variantCount-1
            # Note the STRICT `>`. A variant that merely draws with the
            # incumbent does not displace it. That is not a detail — it is
            # what makes "consistent" mean anything: a descriptor moves only
            # when something is genuinely better, so a scenario in which every
            # descriptor is merely tied for best is still consistent. And
            # among several equally good challengers the lowest-numbered wins,
            # because it is reached first and later equals cannot unseat it.
            if descriptorScores[variantIndex + 1] > bestScore
                bestScore = descriptorScores[variantIndex + 1]
                successor[descriptorIndex] = variantIndex
            end
        end

        firstColumn = lastColumn + 1   # on to the next descriptor's block
    end
    return successor
end

function succession_step(::SequentialSuccession, cib::CIB, scenario::Vector{Int})
    successor = copy(scenario)

    @inbounds for descriptorIndex in 1:cib.numberOfDescriptors
        # THE ONE DIFFERENCE from the global rule: the scores are recomputed
        # here, inside the loop, from `successor` — the scenario as it has
        # already been changed by the descriptors that moved earlier in this
        # same step. So if descriptor 1 has just switched variants, descriptor
        # 2 now scores against a world in which that switch has happened.
        #
        # Two consequences. The rule is sensitive to the order the descriptors
        # happen to be listed in, and a step costs one full rescoring per
        # descriptor instead of one per step.
        impactBalance = impact_balance(cib, successor)

        # The same contest as the global rule, but reaching into the long
        # score vector directly rather than walking a sliding window:
        # descriptor i's block starts at desc_offsets[i], so its variant j
        # sits at desc_offsets[i] + j + 1. Two spellings of one arithmetic.
        offset = cib.desc_offsets[descriptorIndex]
        variantCount = cib.numberOfVariants[descriptorIndex]

        # The incumbent seeds the contest again, so ties keep it in place.
        bestScore = impactBalance[offset + successor[descriptorIndex] + 1]
        for variantIndex in 0:variantCount-1
            if impactBalance[offset + variantIndex + 1] > bestScore
                bestScore = impactBalance[offset + variantIndex + 1]
                successor[descriptorIndex] = variantIndex
            end
        end
    end
    return successor
end

# Convenience method: calling succession_step without naming a rule means the
# standard global rule. (This one-line `f(x) = ...` form is just compact
# syntax for defining a function — the body is the expression after the `=`.)
succession_step(cib::CIB, scenario::Vector{Int}) =
    succession_step(GlobalSuccession(), cib, scenario)

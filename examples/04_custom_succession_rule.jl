"""
Plugging in a custom succession dynamics.

The succession *rule* — the deterministic map that sends a scenario to its
successor — is an extension point. To add a new one you subtype
`SuccessionRule` and define a single method, `succession_step(rule, cib, u)`.
Every analysis routine (`succession`, `find_consistent`, `find_basins`) then
works with it immediately through a generic, rule-agnostic path; the built-in
`GlobalSuccession` additionally carries the fast threaded implementations.

Run from the repo root:
    julia --project=. examples/04_custom_succession_rule.jl
"""

using CrossImpactBalances

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files", "CIB_global.scw")
cib = load_scw(SAMPLE)

# ── A custom rule: "inertial" global succession ──────────────────────────────
# Like global succession, but a descriptor only moves off its current variant
# when some alternative beats it by a strict margin `m`. Recovers ordinary
# global succession at m = 0. This is a toy example; the point is that the
# whole file below is all the code a new dynamics needs.
struct InertialSuccession <: SuccessionRule
    margin::Int
end

function CrossImpactBalances.succession_step(rule::InertialSuccession,
                                             cib::CIB, u::Vector{Int})
    ib = impact_balance(cib, u)
    v = copy(u)
    for i in 1:cib.ndesc
        off = cib.desc_offsets[i]
        nv = cib.nvariants[i]
        cur = ib[off + u[i] + 1]
        best = u[i]
        bestval = cur
        for j in 0:nv-1
            # require a strict margin over the *current* variant's score
            if ib[off + j + 1] > bestval && ib[off + j + 1] - cur > rule.margin
                bestval = ib[off + j + 1]
                best = j
            end
        end
        v[i] = best
    end
    return v
end

# ── Use it exactly like the built-in rule ────────────────────────────────────
for rule in (GlobalSuccession(), InertialSuccession(1), InertialSuccession(3))
    kern = find_consistent(cib; rule=rule)
    _, sizes, cycles = find_basins(cib; rule=rule)
    println(rpad(string(rule), 24),
            "consistent scenarios: ", length(kern),
            "   | largest basin: ", isempty(sizes) ? 0 : maximum(sizes),
            "   | into-cycle starts: ", cycles)
end

# A larger margin makes more scenarios "sticky" and therefore consistent,
# which is exactly the kind of behavioural change a new succession algorithm
# can express — with no modification to the search or basin engines.

# ── Levers: which commitment moves one future into another ──────────────────
#
# Basin sizes say how much of the possibility space drains into each future.
# This asks the question a decision-maker actually has: we are HERE, and we
# can change one thing — which futures does that open, and which does it shut?
#
# The model below is deliberately small enough to check by hand. Five
# descriptors, two variants each, wired so that three things are true:
#
#   * Trade and Economy reinforce each other (free trade → growth → free
#     trade), so that pair is bistable;
#   * Policy and Emissions do the same, but with a twist that lets them chase
#     each other in a loop instead of settling;
#   * Worldview votes on Trade but nobody votes on Worldview — an exogenous
#     dial, and the interesting question is how much it decides.
#
# Nothing here scales with the scenario space: every walk starts from a known
# attractor. The same code runs unchanged on a model with 10^18 scenarios,
# which is exactly why this analysis exists.

using CrossImpactBalances

descriptors = ["Worldview", "Trade", "Economy", "Policy", "Emissions"]
variants = Dict("Worldview" => ["Cooperative", "Fragmented"],
                "Trade"     => ["Free", "Protectionist"],
                "Economy"   => ["Growth", "Stagnation"],
                "Policy"    => ["Ambitious", "Weak"],
                "Emissions" => ["Falling", "Rising"])
variantCounts = [2, 2, 2, 2, 2]
offsets = cumsum(vcat(0, variantCounts[1:end-1]))
cim = zeros(Int, 10, 10)

# Worldview → Trade. Nobody votes back, so this is the dial.
cim[1, 3], cim[1, 4] = 3, -3        # Cooperative pushes Free trade
cim[2, 3], cim[2, 4] = -3, 3        # Fragmented pushes Protectionism
# Trade ↔ Economy, mutually reinforcing.
cim[3, 5], cim[3, 6] = 2, -2        # Free trade → Growth
cim[4, 5], cim[4, 6] = -2, 2        # Protectionism → Stagnation
cim[5, 3], cim[5, 4] = 1, -1        # Growth → Free trade
cim[6, 3], cim[6, 4] = -1, 1        # Stagnation → Protectionism
# Policy ↔ Emissions, also mutually reinforcing. Note what this produces
# under the standard rule, where every descriptor moves at once: the two
# ALIGNED combinations settle, and the two anti-aligned ones (ambitious policy
# with rising emissions, weak policy with falling emissions) flip both
# descriptors simultaneously forever — a cycle, and a good illustration that
# a pair can be perfectly sensible and still never come to rest.
cim[7, 9], cim[7, 10] = 2, -2       # Ambitious → Falling
cim[8, 9], cim[8, 10] = -2, 2       # Weak → Rising
cim[9, 7], cim[9, 8] = 1, -1        # Falling sustains Ambitious
cim[10, 7], cim[10, 8] = -1, 1      # Rising sustains Weak

cib = CIB(descriptors, variants, variantCounts, cim, permutedims(cim), 10, 5,
          Vector{Vector{Int}}(), offsets)

println("Scenario space: ", scenario_count(cib), " scenarios")
println("Consistent scenarios: ", length(find_consistent(cib)))

# ── The lever map between futures ───────────────────────────────────────────
println("\n", "="^70)
graph = transition_graph(cib)
show(stdout, MIME"text/plain"(), graph)
println()

# ── Add the current state of the world ──────────────────────────────────────
# A starting point that need not be consistent — where things stand today.
# Its first edge is the baseline: where we drift if nobody does anything.
println("\n", "="^70)
today = Dict("Worldview" => "Fragmented",
             "Trade"     => "Free",
             "Economy"   => "Growth",
             "Policy"    => "Weak",
             "Emissions" => "Rising")
withWorld = transition_graph(cib; from=today)
show(stdout, MIME"text/plain"(), withWorld)
println()

# ── Combination levers ──────────────────────────────────────────────────────
# radius=2 also tries pairs of changes, and reports only those reaching
# somewhere no single change could — the futures you cannot get to one step
# at a time.
println("\n", "="^70)
paired = transition_graph(cib; from=today, radius=2)
pairLevers = [edge for edge in paired.edges if length(edge.changes) == 2]
println("Combination levers found (unreachable by any single change): ", length(pairLevers))
for edge in pairLevers
    changeText = join([string(paired.descriptors[d], "=",
                              paired.variantNames[d][v + 1]) for (d, v) in edge.changes], " + ")
    destination = join([paired.variantNames[i][paired.nodes[edge.to][i] + 1]
                        for i in 1:5], "/")
    println("  ", changeText, "  ⟶  ", destination,
            paired.kinds[edge.to] === :cycle ? "  (a cycle)" : "")
end

# ── Render it ───────────────────────────────────────────────────────────────
# Paste into any Graphviz viewer, or:  dot -Tsvg levers.dot -o levers.svg
println("\n", "="^70)
println("Graphviz source (first lines):")
for line in first(split(to_dot(withWorld), '\n'), 12)
    println("  ", line)
end
println("  …")

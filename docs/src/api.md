# API reference

## Overview

```@docs
CrossImpactBalances
```

## Loading data

```@docs
load_scw
load_solutions
CIB
```

## Scenario encoding

```@docs
signature
inv_signature
max_signature
scenario_count
```

## Scoring

```@docs
impact_balance
```

## Editing a loaded model

Both matrices are edited in place, so a model can be re-analysed without
re-parsing its `.scw` file — which is what makes sensitivity sweeps cheap.

```@docs
set_impact!
get_impact
```

## Succession rules

The succession dynamics is pluggable: define a `struct MyRule <: SuccessionRule`
and a single `succession_step(::MyRule, cib, u)` method, and every analysis
routine works with it through a generic path.

```@docs
SuccessionRule
GlobalSuccession
SequentialSuccession
CrossImpactBalances.fixed_point_margin
```

## Dynamics

```@docs
succession_step
find_consistent
find_basins
```

## Very large scenario spaces

For spaces beyond the table method's memory (see the README's "Very large
scenario spaces" for the decision table): exact streaming lives on
`find_basins` itself (`method=:stream`, `signature_range`, `cache_bytes`);
share estimation and structural decomposition are their own entry points.

```@docs
estimate_basins
BasinEstimate
influence_structure
InfluenceStructure
fix_descriptor
split_cib
product_basins
ComposedBasins
```

## Levers between futures

Basin sizes measure volume; this measures connectivity — which commitment
moves the system from one consistent scenario into another, and where the
current state of the world drifts if nothing is done. Every walk starts from a
known attractor, so the cost does not depend on the size of the scenario space
at all.

```@docs
transition_graph
TransitionGraph
Transition
to_dot
```

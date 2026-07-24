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
```

## Scoring

```@docs
impact_balance
```

## Succession rules

The succession dynamics is pluggable: define a `struct MyRule <: SuccessionRule`
and a single `succession_step(::MyRule, cib, u)` method, and every analysis
routine works with it through a generic path.

```@docs
SuccessionRule
GlobalSuccession
SequentialSuccession
```

## Dynamics

```@docs
succession_step
find_consistent
find_basins
```

## Threshold-gated fluctuation analysis

```@docs
sim_anneal
build_graph
merge_scenarios
```

## In-place editing

Tweak an already-loaded model's cross-impact matrix without re-parsing the
`.scw` file. Both the matrix and its stored transpose are updated together, so
every analysis routine stays correct.

```@docs
set_impact!
get_impact
```

## Utilities

```@docs
set_thresholds!
rand_scenario
CIB
```

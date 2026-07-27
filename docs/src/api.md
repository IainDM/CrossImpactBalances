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
CrossImpactBalances.fixed_point_margin
```

## Dynamics

```@docs
succession_step
find_consistent
find_basins
```

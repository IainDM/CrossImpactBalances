# CrossImpactBalances.jl

Julia implementation of Cross-Impact Balance (CIB) scenario analysis — a
methodology for finding the internally consistent futures of a qualitative
multi-descriptor system from an expert-elicited cross-impact matrix.

This package is a reimplementation of the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with significantly improved performance: the threaded exhaustive search and
memoized basin-of-attraction analysis routinely deliver 100× speedups on the
benchmark suite (see `test/benchmark.jl`).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/IainDM/CrossImpactBalances.jl")
```

For multi-threaded exhaustive search, start Julia with multiple threads:

```bash
julia -t auto --project=.
```

## Quick start

```julia
using CrossImpactBalances

cib = load_scw("test/sample_files/CIB_global.scw")

# Find all consistent scenarios (fixed points of the succession operator)
for u in cib.kernel
    println("Scenario signature ", signature(cib, u))
    for (i, desc) in enumerate(cib.descriptors)
        println("  ", desc, " = ", cib.variants[desc][u[i] + 1])
    end
end

# Exhaustive basin-of-attraction analysis
fixed_points, basin_sizes, cycle_count = find_basins(cib)
```

For larger problems, use `load_scw(path; exhaustive=true)` to enumerate the
full scenario space in parallel across threads.

## What's included

| Function | Purpose |
|----------|---------|
| `load_scw`, `load_solutions` | Parse ScenarioWizard `.scw` and `.sl` files |
| `signature`, `inv_signature`, `max_signature` | Bijection between scenarios and integers |
| `impact_balance`, `own_impact_balance`, `cross_impact_balance`, `inner_product` | CIB scoring primitives |
| `succession_step`, `succession` | Deterministic global-succession dynamics |
| `find_consistent` | Find all fixed points (Monte-Carlo or exhaustive) |
| `find_basins` | Exhaustive basin-of-attraction analysis with memoization (Julia-only addition) |
| `sim_anneal`, `build_graph`, `merge_scenarios` | Threshold-gated fluctuation analysis for kernel reduction |
| `inner_product_matrix` | Pairwise similarity of kernel scenarios |

## Citation

If you use this software, please cite the CIB methodology:

> Weimer-Jehle, W. (2006). *Cross-impact balances: A system-theoretical
> approach to cross-impact analysis*. Technological Forecasting and Social
> Change, 73(4), 334–361.

## Acknowledgments

The algorithms in this package are a port of the Python implementation by the
Stockholm Environment Institute,
[sei-international/cibsa](https://github.com/sei-international/cibsa), with
performance optimizations and a new basin-analysis routine.

## License

MIT — see [LICENSE](LICENSE). Original CIBSA copyright is preserved.

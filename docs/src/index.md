# CrossImpactBalances.jl

Julia implementation of Cross-Impact Balance (CIB) scenario analysis — a
methodology for finding the internally consistent futures of a qualitative
multi-descriptor system from an expert-elicited cross-impact matrix.

This package is a reimplementation of the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with a memoized basin-of-attraction analysis and a threaded exhaustive search.
The benchmark suite (`test/benchmark.jl`) routinely reports 25–100× speedups
over the Python reference.

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

cib = load_scw("path/to/CIB_global.scw")

# Each consistent scenario
for u in cib.kernel
    println("sig=", signature(cib, u))
    for (i, desc) in enumerate(cib.descriptors)
        println("  ", desc, " = ", cib.variants[desc][u[i] + 1])
    end
end

# Exhaustive basin-of-attraction analysis
fixed_points, basin_sizes, cycle_count = find_basins(cib)
```

For more, see the [`examples/`](https://github.com/IainDM/CrossImpactBalances.jl/tree/main/examples)
folder in the repository and the [API reference](@ref).

## Citation

```bibtex
@article{weimer-jehle-2006-cib,
  author  = {Weimer-Jehle, Wolfgang},
  title   = {Cross-impact balances: A system-theoretical approach to cross-impact analysis},
  journal = {Technological Forecasting and Social Change},
  volume  = {73},
  number  = {4},
  pages   = {334--361},
  year    = {2006},
}
```

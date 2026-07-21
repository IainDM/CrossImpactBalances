# CrossImpactBalances.jl

Julia implementation of Cross-Impact Balance (CIB) scenario analysis — a
methodology for finding the internally consistent futures of a qualitative
multi-descriptor system from an expert-elicited cross-impact matrix.

This package is a reimplementation of the Python
[sei-international/cibsa](https://github.com/sei-international/cibsa) library
with substantial performance work and one new analysis routine. See
[Performance](#performance) for measured speedups across problem sizes.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/IainDM/CrossImpactBalances.jl")
```

For multi-threaded exhaustive search and basin analysis, start Julia with
multiple threads:

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
| `find_consistent` | Find all fixed points (Monte-Carlo, or `exhaustive=true` with `algorithm=:auto`/`:bnb`/`:sweep`) |
| `find_basins` | Two-phase basin-of-attraction analysis (Julia-only addition) |
| `sim_anneal`, `build_graph`, `merge_scenarios` | Threshold-gated fluctuation analysis for kernel reduction |
| `inner_product_matrix` | Pairwise similarity of kernel scenarios |
| `set_thresholds!`, `rand_scenario` | API helpers matching the Python reference |

## Performance

Same machine (4-core Intel Cascade Lake), Julia 1.12 with `-t auto`,
Python 3 single-process. Median of three runs.

### Same algorithm — language + SIMD speedup

| Benchmark | Scenarios | Python `find_consistent` | Julia `find_consistent` | Speedup |
|---|---:|---:|---:|---:|
| `bench_medium` | 1,024 | 71 ms | 1.6 ms | **44×** |
| `bench_large` | 6,561 | 843 ms | 37 ms | **23×** |
| `bench_xlarge` (MC sample, 10k of 59k) | 59,049 | 1.92 s | 47 ms | **41×** |
| `bench_xlarge` (full enumeration) | 59,049 | 9.79 s | 47 ms | **210×** |
| `bench_50x50` (MC sample) | 60,466,176 | 21.7 s | 371 ms | **59×** |

### Julia-only fast paths

The specialized exhaustive routine (threaded + mixed-radix counter +
per-descriptor early exit) and the memoized basin analysis are new in
`CrossImpactBalances.jl`. They have no Python equivalent — the Python
column is what you'd get from `find_consistent` (raised threshold) and a
naive Python basin loop.

| Operation | Python (naive equivalent) | Julia | Speedup |
|---|---:|---:|---:|
| `bench_xlarge` full enumeration (`exhaustive=true`) | 9.8 s | **1.3 ms** | **~7,500×** |
| `bench_xlarge` basin analysis (`find_basins`) | 10.0 s | **8.1 ms** | **~1,200×** |
| `bench_typical` (59k scenarios, same problem class as xlarge) | ~10 s | **1.3 ms** | **~7,700×** |

The basin numbers stack four wins multiplicatively: language (~30×) ×
memoization (~10×) × threading (~2.3× on 4 cores) × row-major SIMD (~2×).

Since those measurements, `find_basins` was rebuilt around a single label
array shared by all threads (4 bytes per scenario instead of 8 bytes per
scenario *per thread*), narrow-integer SIMD scoring, work-stealing block
scheduling, and several interleaved walks per thread to overlap the random
label-array reads that otherwise bound the walk on DRAM latency. On a 4-core
AVX-512 box this makes basin analysis another ~2× faster on the 60M-scenario
model and ~6× faster on the 408M-scenario `CIB_nested` model — which
previously could not run multi-threaded at all within 16 GB of RAM (the old
per-thread caches would have needed ~13 GB; the shared array needs 1.6 GB).

To reproduce: `julia --project=. -t auto test/benchmark.jl`,
`test/benchmark_50x50.jl`, and `test/benchmark_basins.jl nested`.

### v0.2 optimizations

Version 0.2 rewrites all three engines (results are bit-identical; a
property-test suite pins every algorithm against brute-force oracles):

- **Exhaustive sweep**: the impact balance is maintained *incrementally* as a
  mixed-radix odometer walks the space (two-row `@simd` delta per step, Int16
  accumulation when the matrix allows), replacing per-scenario score
  recomputation; work is split into fine-grained tasks for load balancing.
- **`find_basins`**: two phases — a threaded flat successor table, then a
  path-compressed resolution walk with O(1) cycle detection. Memory drops
  from `nthreads × 8n` bytes to a thread-count-independent `~8n` bytes.
- **Branch-and-bound** (`algorithm=:bnb`, chosen automatically by
  `algorithm=:auto` for spaces ≥ 10^5): descriptors are assigned depth-first
  and any subtree in which some assigned descriptor's chosen variant is
  provably beaten in every completion is discarded using precomputed suffix
  score bounds. Exact, and typically visits a few percent of the space; a
  visited-node budget falls back to the sweep on weakly-coupled matrices
  where pruning cannot pay off.

Measured on a 4-core Intel Xeon @ 2.10 GHz, Julia 1.11.7 `-t 4`, find-only
medians of 3 (`test/bench_optim.jl`):

| File | Scenarios | v0.1 sweep | v0.2 sweep | v0.2 B&B (nodes visited) | v0.1 basins | v0.2 basins |
|---|---:|---:|---:|---:|---:|---:|
| `bench_typical` | 59,049 | 0.59 ms | 0.39 ms | 0.22 ms (13.3%) | 5.8 ms | 2.2 ms |
| `bench_xlarge` | 59,049 | 0.61 ms | 0.39 ms | 0.22 ms (15.2%) | 5.5 ms | 1.8 ms |
| `bench_50x50` | 60,466,176 | 915 ms | 348 ms | **23.6 ms** (2.25%) | 11.0 s | **1.91 s** |

Net effect on the 60M-scenario stress file: exhaustive search is ~39× faster
than v0.1 (and B&B's cost scales with the pruned tree, not the space, so the
gap widens on larger problems); full basin analysis is ~5.8× faster.

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

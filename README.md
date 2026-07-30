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
for u in cib.consistentScenarios
    println("Scenario signature ", signature(cib, u))
    for (i, desc) in enumerate(cib.descriptors)
        println("  ", desc, " = ", cib.variants[desc][u[i] + 1])
    end
end

# Exhaustive basin-of-attraction analysis
fixed_points, basin_sizes, cycle_count = find_basins(cib)
```

`load_scw` always enumerates the full scenario space in parallel across
threads — there is no sampling mode and no completeness caveat. Pass
`algorithm=:sweep` or `:bnb` to force a search strategy; the default `:auto`
picks one from the space size.

The package has no dependencies outside the Julia standard library — in fact,
no dependencies at all.

## Python interface

A Pythonic wrapper lives in [`python/`](python/README.md). It embeds this Julia
engine **in-process** via [`juliacall`](https://github.com/JuliaPy/PythonCall.jl),
so a parsed model stays resident — you can batch-run many `.scw` files and even
tweak individual cross-impact values in place without re-parsing.

```bash
pip install -e "python/[pandas]"          # from the repo root
```

```python
from crossimpactbalances import Model, run_models, sweep_impact

m = Model.load("test/sample_files/CIB_global.scw")
m.consistent_scenarios()                  # [{'WTRD': 'FT', 'WSEC': 'Rlx', ...}, ...]

# Edit an expert judgement and re-run — no re-parse of the .scw.
m.set_impact(source=("WTRD", "FT"), target=("WSEC", "Alrt"), value=3)
m.consistent_scenarios()

# Or sweep one impact / batch over many model files.
sweep_impact(m, ("WTRD", "FT"), ("WSEC", "Alrt"), [-3, 0, 3, 6])
run_models("test/sample_files", analysis="basins")
```

Set `PYTHON_JULIACALL_THREADS=auto` before first use for multi-threaded
exhaustive/basin analysis. The two new engine primitives that back in-place
editing — `set_impact!` and `get_impact` — are also exported for Julia callers.
See [`python/README.md`](python/README.md) and
[`examples/python/run_many.py`](examples/python/run_many.py).

### Two backends (including a source-free one)

The same `Model` API runs on either backend, chosen by `Model.load(...,
backend=...)` or the `CIB_BACKEND` env var:

- **`juliacall`** — runs the Julia source package in-process (above). Best for
  developing in this repo.
- **`native`** — a compiled `libcib` shared library (machine code + bundled
  Julia runtime) driven via `ctypes`, built with PackageCompiler from
  [`capi/`](capi/README.md). It ships **no Julia source** and needs no Julia
  install, so it's the way to distribute the engine as a black box. Build it
  with `julia --project=build build/build_library.jl`; see
  [`python/README.md`](python/README.md) for bundling it into a wheel.

## Desktop app

A point-and-click version for non-programmers lives in [`app/`](app): browse to
a `.scw` file, then **Find Consistent Scenarios** or **Find Basins**, view the
results (with basin sizes) and export the full basin analysis as CSV. It runs
locally in your browser and has no external dependencies. Run it from source
with `julia --project=app -t auto app/run.jl`, or build a standalone Windows
installer (`.msi`) — no Julia needed on the target machine — via the pipeline
in [`build/`](build/README.md).

## What's included

| Function | Purpose |
|----------|---------|
| `load_scw`, `load_solutions` | Parse ScenarioWizard `.scw` and `.sl` files |
| `signature`, `inv_signature`, `max_signature` | Bijection between scenarios and integers |
| `impact_balance` | Score every variant against a scenario |
| `set_impact!`, `get_impact` | Read/edit one cross-impact in a loaded model, no re-parse |
| `SuccessionRule`, `GlobalSuccession`, `SequentialSuccession` | Pluggable succession dynamics (see [Pluggable succession rules](#pluggable-succession-rules)) |
| `succession_step` | One deterministic succession step under a chosen rule |
| `find_consistent` | Find every fixed point by exhaustive search (`algorithm=:auto`/`:bnb`/`:sweep`); honours `rule=` |
| `find_basins` | Two-phase basin-of-attraction analysis (Julia-only addition); honours `rule=` |

## Pluggable succession rules

The succession *dynamics* — the deterministic map from a scenario to its
successor — is an extension point. `GlobalSuccession` (the classical
ScenarioWizard/CIBSA rule) is the default and carries the fast threaded
sweep, branch-and-bound, and two-phase basin implementations. To add a new
rule, subtype `SuccessionRule` and define one method:

```julia
struct MyRule <: SuccessionRule end

function CrossImpactBalances.succession_step(::MyRule, cib::CIB, u::Vector{Int})
    # return the successor scenario (0-based variant indices)
end
```

That's the entire contract. Every analysis routine then works with it:

```julia
find_consistent(cib; rule=MyRule())
find_basins(cib; rule=MyRule())
succession_step(MyRule(), cib, u)
```

A custom rule runs through a generic path that calls `succession_step` once
per scenario (a correctness baseline — threaded, but slower than the
odometer fast path the default rule gets). A rule whose fixed points can be
expressed as "no variant beats the current one by more than a margin" can
declare `CrossImpactBalances.fixed_point_margin(rule)` and inherit the fast
sweep and branch-and-bound searches for free.
`SequentialSuccession` (Gauss–Seidel-style updates)
ships as a second built-in rule, and
[`examples/04_custom_succession_rule.jl`](examples/04_custom_succession_rule.jl)
is a worked example. The `algorithm=:sweep/:bnb` search-strategy switch
applies only to rules with a `fixed_point_margin` — both built-in rules
declare one (`GlobalSuccession` trivially; `SequentialSuccession` because
its fixed points provably coincide with global's, even though its
trajectories differ).

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
  provably beaten in every completion is discarded using precomputed
  per-variant-pair suffix difference bounds — for each (chosen, rival) pair
  the extreme of the score *difference* is taken over each undecided
  descriptor's shared variant choice, which is strictly tighter than
  bounding the two scores independently. Exact, and typically visits a few
  percent of the space; a visited-node budget falls back to the sweep on
  weakly-coupled matrices where pruning cannot pay off.

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

A subsequent tightening of the branch-and-bound test — bounding each
(chosen, rival) score *difference* over the undecided descriptors' shared
variant choices, instead of bounding the two scores independently — shrinks
the visited tree further, with bit-identical results: `bench_typical`
13.3% → 6.66% of the space, `bench_xlarge` 15.2% → 8.71%, `bench_50x50`
2.25% → 0.601% (node counts are machine-independent; the timings in the
table above predate this change). On `bench_50x50` the smaller tree makes
the B&B search ~4× faster again.

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

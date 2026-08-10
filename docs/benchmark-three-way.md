# Three-way CIB benchmark: ScenarioWizard vs CIBSA vs JuCIB

## Goal

Measure wall-clock time and verify correctness for three implementations of
Cross-Impact Balance scenario analysis on a common set of `.scw` files:

1. **ScenarioWizard** — Weimer-Jehle's reference desktop app (Windows).
2. **CIBSA** — Python implementation, [sei-international/cibsa](https://github.com/sei-international/cibsa).
3. **JuCIB** — `CrossImpactBalances.jl`, [IainDM/CrossImpactBalances](https://github.com/IainDM/CrossImpactBalances) (branch `claude/cim-rowmajor`).

The headline metric is **time to find all consistent scenarios** ("kernel"
size and contents). Where the tool supports it, also report time to compute
**basin sizes** for every scenario.

## Hardware / config to log

Before benchmarking, capture and include in the report:

- CPU model, physical core count, total RAM.
- OS and version.
- ScenarioWizard version (Help → About).
- Python version and `numpy` version.
- Julia version (`julia --version`).
- Julia thread count (`Threads.nthreads()` in REPL, or the value of
  `JULIA_NUM_THREADS` / `julia -t auto`).

## Test files

Run all three implementations on the **same** `.scw` files. Suggested set,
small to large:

| File | ~ Scenarios | Where to get it |
|---|---:|---|
| `CIB_global.scw` | 36 | `test/sample_files/` in the JuCIB repo |
| `bench_typical.scw` | 59,049 | same |
| `bench_xlarge.scw` | 59,049 | same |
| `bench_50x50.scw` | 60,466,176 | same |
| the Weimer-Jehle corpus | 10⁸·⁴–10³⁰·¹ | 38 models from the method's author; not redistributed — see [`test/WWJ_CORPUS.md`](../test/WWJ_CORPUS.md) |

The small files are sanity checks (all three implementations should agree
on the kernel). The mid-sized files measure same-algorithm speedup. The 60 M
file stress-tests the big-problem regime. The Weimer-Jehle corpus is the
correctness argument rather than the speed one: 19 of its models ship with
ScenarioWizard's own solution set, at sizes no other tool here can reach.

## Running each implementation

### ScenarioWizard

1. Launch ScenarioWizard, open the `.scw` file (File → Open).
2. Find the menu/button that runs the consistency search. (Likely under
   *Analysis* or similar — exact label depends on version.)
3. **Time it manually with a stopwatch.** Start the clock when you click
   "search"; stop when the kernel/results dialog appears.
4. Record three trials; report the median.
5. Record the number of consistent scenarios found.
6. If ScenarioWizard reports its own elapsed time anywhere, capture that too.

Notes:
- ScenarioWizard is GUI-driven, so 100 ms timing resolution isn't possible.
  Round to the nearest 0.5 s for small problems and the nearest second for
  large ones.
- If the search has a configurable algorithm (e.g. exhaustive vs.
  sampled), use the **default** for a first comparison and note the
  setting in the report.

### CIBSA (Python)

```bash
git clone https://github.com/sei-international/cibsa
cd cibsa
pip install numpy
```

Time `find_consistent` with `mc_threshold` set high enough to force full
enumeration on the smaller files. For files larger than the default
`mc_threshold=10000`, run **once** with MC (matching ScenarioWizard's
default behaviour if any) and **once** with `mc_threshold=10**9` to force
full enumeration. Report both.

```python
import time
from CIB_sim_anneal import CIB

# Full enumeration
cib = CIB("/path/to/file.scw", mc_threshold=10**9)
trials = []
for _ in range(3):
    t0 = time.perf_counter()
    cib.find_consistent()
    trials.append(time.perf_counter() - t0)
print(f"median: {sorted(trials)[1]:.3f}s  kernel: {len(cib.kernel)}")
```

For the 4 M Cascade file the full-enumeration Python run may take many
hours — if so, abandon it after ~30 minutes and report "did not finish"
rather than waiting.

### JuCIB (Julia)

```bash
git clone -b claude/cim-rowmajor https://github.com/IainDM/CrossImpactBalances
cd CrossImpactBalances
julia -t auto --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'   # one-off
```

Time *three* calls per file: the sweep (the like-for-like against CIBSA's full
enumeration — every scenario is checked), branch-and-bound, and `find_basins`.
Run with **the same thread count as the host has physical cores** unless
investigating scaling.

```julia
using CrossImpactBalances
scw = "/path/to/file.scw"
cib = load_scw(scw; kernel=Vector{Vector{Int}}())   # parse only, no search

# Warm-up (eliminates JIT compile time)
find_basins(cib)

# Check every scenario — the like-for-like against CIBSA's full enumeration.
# Only available while the model has fewer than typemax(Int64) scenarios.
ts = Float64[]
for _ in 1:3
    t0 = time_ns()
    find_consistent(cib; algorithm=:sweep)
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_sweep: $(ts[2]) s")

# Branch-and-bound: identical answer, most of the space never visited
ts = Float64[]
for _ in 1:3
    t0 = time_ns()
    find_consistent(cib; algorithm=:bnb)
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_bnb:   $(ts[2]) s")

# Basin analysis
ts = Float64[]
for _ in 1:3
    t0 = time_ns()
    find_basins(cib)
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_bas:   $(ts[2]) s")
```

`load_scw(scw)` on its own also runs `find_consistent` (with `algorithm=:auto`),
which is what an ordinary caller does; pass `kernel=Vector{Vector{Int}}()` as
above when you want to time the search separately from the parse.

Use `julia -t auto` (or `-tN` for N threads) and **set `JULIA_NUM_THREADS`
to match the host's physical core count** for the headline runs.

## Verification

Before claiming any speedup is real, verify the three implementations
**agree on the kernel**. For each `.scw`:

- Print the list of consistent scenarios from each tool (descriptor →
  variant choices for each).
- Convert to a comparable form. **Sorted tuples of variant indices, not
  signatures**: `signature` accumulates in `Int` and wraps silently once a model
  passes `typemax(Int64)` scenarios, so two different scenarios can share a
  value. `sort` on `Vector{Vector{Int}}` is lexicographic and total, so
  `sort(a) == sort(b)` is exact set equality at any model size.
- Report the two directions **separately**. "Present in the reference but not
  found" means the search is unsound and is the most serious result available;
  "found but not in the reference" is either a bug or a finding about the other
  tool. A single agree/disagree boolean hides the difference.
- Note any disagreements. Possible legitimate reasons:
  - CIBSA may sample (MC) and miss small-basin fixed points; JuCIB always
    enumerates the whole space and will find them all.
  - Different tie-breaking rules in `succession_step` — *should* be the
    same (favour current variant, then lowest index), but worth checking
    on small problems.
  - The reference `.sl` was exported from an earlier revision of the matrix.
    This is not hypothetical: two files in the Weimer-Jehle corpus are exactly
    that, and load cleanly while agreeing on nothing.

If there are disagreements that aren't explained by sampling, that's a
finding, not a measurement bug. Flag it.

### External reference corpus

`test/verify_wwj.jl` runs the whole protocol above against Weimer-Jehle's 38
models and writes `test/bench_results_wwj.json`, exiting non-zero on any
discrepancy. It needs the corpus staged at `test/wwj_corpus/` — see
[`test/WWJ_CORPUS.md`](../test/WWJ_CORPUS.md). Results are summarised in
[`BENCHMARK_THREEWAY_RESULTS.md`](../BENCHMARK_THREEWAY_RESULTS.md).

## Reporting

Produce one table per benchmark file:

| Implementation | Mode | Median time | Kernel size | Notes |
|---|---|---:|---:|---|
| ScenarioWizard | default | _ s | _ | stopwatched |
| CIBSA | full enum | _ s | _ | |
| CIBSA | MC sample | _ s | _ | mc_threshold=10000 |
| JuCIB | sweep (every scenario) | _ s | _ | `algorithm=:sweep`, -t auto |
| JuCIB | branch-and-bound | _ s | _ | `algorithm=:bnb`, -t auto |
| JuCIB | find_basins | _ s | n/a | -t auto |

Then a single summary table across files with the headline speedups,
plus a paragraph noting:

- Whether the kernels agreed across implementations.
- Whether ScenarioWizard's "much faster than CIBSA" reputation holds up
  against JuCIB.
- Any file ScenarioWizard couldn't handle (memory, time-out, etc.).

## Three implementations, one repo

If it helps, drop the raw timings into a JSON file alongside the report so
they can be plotted or rolled into JuCIB's `test/benchmark.jl` later.

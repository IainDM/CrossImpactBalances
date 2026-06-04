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
| (Cascade Institute file) | ~4 M | bring your own |

The small files are sanity checks (all three implementations should agree
on the kernel). The mid-sized files measure same-algorithm speedup. The 60 M
and 4 M files stress-test the big-problem regime.

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

Time *three* calls per file: same-algorithm `find_consistent` (matches
CIBSA), specialised `find_consistent(..., exhaustive=true)`, and
`find_basins`. Run with **the same thread count as the host has physical
cores** unless investigating scaling.

```julia
using CrossImpactBalances
scw = "/path/to/file.scw"

# Warm-up (eliminates JIT compile time)
let cib = load_scw(scw); find_basins(cib); end

# Same algorithm as CIBSA (walks succession from each scenario)
ts = Float64[]
for _ in 1:3
    t0 = time_ns()
    load_scw(scw; mc_threshold=10^9)   # force full enumeration
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_find: $(ts[2]) s")

# Specialised exhaustive (Julia-only fast path)
ts = Float64[]
for _ in 1:3
    t0 = time_ns()
    load_scw(scw; exhaustive=true)
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_exh:  $(ts[2]) s")

# Basin analysis
ts = Float64[]
for _ in 1:3
    cib = load_scw(scw; kernel=Vector{Vector{Int}}())
    t0 = time_ns()
    find_basins(cib)
    push!(ts, (time_ns() - t0) / 1e9)
end
sort!(ts); println("jl_bas:  $(ts[2]) s")
```

Use `julia -t auto` (or `-tN` for N threads) and **set `JULIA_NUM_THREADS`
to match the host's physical core count** for the headline runs.

## Verification

Before claiming any speedup is real, verify the three implementations
**agree on the kernel**. For each `.scw`:

- Print the list of consistent scenarios from each tool (descriptor →
  variant choices for each).
- Convert to a comparable form (sorted signatures, or sorted tuples of
  variant indices).
- Note any disagreements. Possible legitimate reasons:
  - ScenarioWizard and CIBSA may sample (MC) and miss small-basin fixed
    points; JuCIB with `exhaustive=true` will find them all.
  - Different tie-breaking rules in `succession_step` — *should* be the
    same (favour current variant, then lowest index), but worth checking
    on small problems.

If there are disagreements that aren't explained by sampling, that's a
finding, not a measurement bug. Flag it.

## Reporting

Produce one table per benchmark file:

| Implementation | Mode | Median time | Kernel size | Notes |
|---|---|---:|---:|---|
| ScenarioWizard | default | _ s | _ | stopwatched |
| CIBSA | full enum | _ s | _ | |
| CIBSA | MC sample | _ s | _ | mc_threshold=10000 |
| JuCIB | same algo (full enum) | _ s | _ | -t auto |
| JuCIB | exhaustive=true | _ s | _ | -t auto |
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

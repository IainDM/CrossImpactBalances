# Three-way CIB benchmark: ScenarioWizard vs CIBSA vs JuCIB

Wall-clock and correctness comparison of three Cross-Impact Balance (CIB)
consistency-search implementations on a common set of `.scw` files.

- **JuCIB** figures are a fresh run of the **v0.2 engine on 2026-07-21**
  (`test/bench_threeway.jl`, this repo). They supersede the original
  2026-06-04 JuCIB numbers, which were taken on the v0.1 engine before the
  incremental-`Int16` sweep, two-phase basin analysis, and branch-and-bound
  search landed.
- **ScenarioWizard** and **CIBSA** figures are carried over from the
  2026-06-04 run on the **same machine**. Neither tool changed between the
  two dates, so the cross-tool comparison stands; they were not re-timed.

Tools under test:

- **ScenarioWizard (SW)** — Weimer-Jehle's reference desktop app (GUI, driven via UI automation).
- **CIBSA** — published Python package [`sei-international/cibsa`](https://github.com/sei-international/cibsa) (`CIB_sim_anneal.py`).
- **JuCIB** — `CrossImpactBalances.jl` (this repo), v0.2.

## TL;DR

1. **All three tools agree on the kernel** for every file where a full search was
   run — identical consistent-scenario sets, down to the signature. The only
   discrepancy is CIBSA in Monte-Carlo mode missing small-basin fixed points on
   the 60 M-scenario file (expected; see Verification). The v0.2 JuCIB engine
   returns the **same kernels and basin sizes** as v0.1 — only faster.
2. **ScenarioWizard's "much faster than CIBSA" reputation holds** — emphatically.
   SW returns *instantly* on files where CIBSA's pure-Python full enumeration
   takes 15–16 s, and handles the 60 M file in ~5 s where CIBSA cannot enumerate
   it at all.
3. **But JuCIB decisively beats ScenarioWizard.** JuCIB's branch-and-bound
   exhaustive search solves the 60 M file in **~0.013 s (13 ms) and finds all 5
   fixed points** (vs SW's ~5 s for the same 5), and is sub-millisecond on the
   mid-size files. Same algorithm, single-thread, JuCIB is **66–159× faster than
   CIBSA** (pure Julia-vs-Python); with its optimised threaded path it is
   **~1,100–37,000×** faster than CIBSA full-enum on the mid files.
4. JuCIB additionally computes full **basin-of-attraction** sizes — the entire
   60 M-scenario space in **~1.6 s** (two-phase successor table) — which neither
   of the other tools exposes directly.

## Environment

| Item | Value |
|---|---|
| Machine | Lenovo laptop, Windows 11 Home 10.0.26200 |
| CPU | Intel Core i7-13700H — 14 physical cores (6 P + 8 E) / 20 logical |
| RAM | 15.7 GB |
| Julia | 1.12.5, run with `-t 10` (`Threads.nthreads() = 10`) |
| Python | 3.13 (3.13.14 on 2026-07-21; CIBSA figures from 3.13.13 on 2026-06-04), numpy 2.4.1, scipy 1.17.0 |
| ScenarioWizard | 5 (file version 1.0.9522.17909), `C:\Program Files (x86)\ScenarioWizard\` |

> Note on threads: the CPU is hybrid (6 performance + 8 efficiency cores). Julia
> ran with 10 threads, the value configured in `JULIA_NUM_THREADS`. Only JuCIB's
> `exhaustive=true` sweep and `find_basins` table-fill are multi-threaded.

## Method — what was timed

- **Find-only, median of 3** (median of 2 on the 60 M file; warm-up run
  discarded; `time_ns`/`perf_counter`). The model is parsed *once* up front and
  excluded from the timed region, so JuCIB and CIBSA are compared on the search
  alone. Parse cost is negligible everywhere (≤ 0.6 ms even for the 60 M file),
  so full-load ≈ find-only.
- **JuCIB, three modes:** (a) `find_consistent(…; exhaustive=false)` single-thread
  full enumeration — the *same algorithm* as CIBSA (a succession walk from every
  scenario); (b) `find_consistent(…; exhaustive=true)` — the optimised path
  (JuCIB's best): the threaded fixed-point **sweep** for spaces below 10⁵
  scenarios, auto-switching to **branch-and-bound** above it; (c) `find_basins` —
  full basin analysis (threaded successor-table fill + sequential path-compressed
  resolve).
- **CIBSA, two modes:** full enumeration (`mc_threshold=10**9`) and Monte-Carlo
  sampling (`mc_threshold=10000`, `np.random.seed(999)`).
- **ScenarioWizard:** GUI-only, no headless mode. Driven via UI automation:
  open file → *Analyse ▸ Consistent scenarios* → read count → *Save* solutions to
  `<name>_sw.sl`. Timing is **coarse** (screenshot-paced wall-clock); SW's whole
  operation (parse + search + render) is one number. For the five small/mid files
  the result appeared instantly with no progress dialog; only the 60 M file showed
  a "Computing…" progress bar.
- **Ground truth = JuCIB exhaustive** (enumerates the whole space → finds every
  fixed point). Every reported scenario from every tool was independently
  re-checked to be a true fixed point (`succession_step(u) == u`).
- The verification bridge uses ScenarioWizard's exported `.sl` files, re-loaded
  into JuCIB (`load_scw(scw; sl_file=…)`) to put SW's kernel in the same signature
  space.

## Results — per file

Times are **seconds, find-only median**. "Kernel" is the number of consistent
scenarios. ✓ = kernel matches the JuCIB-exhaustive reference and all entries
verified as true fixed points. JuCIB rows = v0.2 engine, 2026-07-21;
CIBSA/SW rows = 2026-06-04.

### CIB_global — 36 scenarios (3 descriptors)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.000014 | 4 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000087 | 4 | ✓ |
| JuCIB | find_basins | 10 | 0.000050 | 4 | ✓ |
| CIBSA | full-enum | 1 | 0.001626 | 4 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 4 | ✓ |

Kernel = `{13, 16, 20, 21}`.

### bench_medium — 1,024 scenarios (5×4)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.001168 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000132 | 2 | ✓ |
| JuCIB | find_basins | 10 | 0.000160 | 2 | ✓ |
| CIBSA | full-enum | 1 | 0.151 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{191, 480}`.

### bench_large — 6,561 scenarios (8×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.016824 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000162 | 2 | ✓ |
| JuCIB | find_basins | 10 | 0.000740 | 2 | ✓ |
| CIBSA | full-enum | 1 | 2.668 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{1973, 2548}`.

### bench_typical — 59,049 scenarios (10×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.229 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000441 | 2 | ✓ |
| JuCIB | find_basins | 10 | 0.00254 | 2 | ✓ |
| CIBSA | full-enum | 1 | 16.35 | 2 | ✓ |
| CIBSA | MC (10 k) | 1 | 2.729 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{13785, 13839}`.

### bench_xlarge — 59,049 scenarios (10×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.232 | 3 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000407 | 3 | ✓ |
| JuCIB | find_basins | 10 | 0.00251 | 3 | ✓ |
| CIBSA | full-enum | 1 | 15.24 | 3 | ✓ |
| CIBSA | MC (10 k) | 1 | 2.686 | 3 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 3 | ✓ |

Kernel = `{15921, 18086, 36255}`.

### bench_50x50 — 60,466,176 scenarios (15 descriptors) — the stress test

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | **exhaustive** (B&B)⁶ | 10 | **0.013** | 5 | ✓ (reference) |
| JuCIB | find_basins | 10 | 1.61 | 5 | ✓ |
| JuCIB | full-enum | 1 | not run² | — | — |
| CIBSA | full-enum | 1 | **DNF**³ | — | — |
| CIBSA | MC (10 k) | 1 | 30.13 | **3** | partial⁴ |
| ScenarioWizard | GUI full | — | **≈ 2–6 s**⁵ | 5 | ✓ |

Kernel (exhaustive) = `{1984298, 8693546, 11540522, 19807290, 31030843}`.

¹ Result appeared instantly with no progress dialog → below the ~0.5 s resolution
of GUI timing.
² Single-thread succession-walk enumeration of 60 M scenarios was skipped; naive
extrapolation from `bench_typical` is on the order of ~4 min. JuCIB's threaded
branch-and-bound path (0.013 s) is the real number.
³ Pure-Python full enumeration of 60 M scenarios is infeasible (hours) — not run,
per the benchmark's "abandon after 30 min" rule.
⁴ CIBSA-MC (seed 999, 10 k samples) found 3 of the 5 fixed points. See Verification.
⁵ Coarse — a "Computing…" progress dialog appeared and the result rendered within
the next ~5 s.
⁶ `exhaustive=true` auto-selects branch-and-bound above 10⁵ scenarios; on this
file it visits **1,361,846 nodes = 2.25 %** of the space. Across runs the median
sat at 11–16 ms; 0.013 s is representative.

## Cross-file summary (find-only, seconds)

| File | Scenarios | Kernel | CIBSA full | CIBSA MC | SW | JuCIB 1-thread | JuCIB exhaustive (10t) | JuCIB basins (10t) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CIB_global | 36 | 4 | 0.0016 | — | instant | 0.000014 | 0.000087 | 0.000050 |
| bench_medium | 1,024 | 2 | 0.151 | — | instant | 0.001168 | 0.000132 | 0.000160 |
| bench_large | 6,561 | 2 | 2.668 | — | instant | 0.016824 | 0.000162 | 0.000740 |
| bench_typical | 59,049 | 2 | 16.35 | 2.73 | instant | 0.229 | 0.000441 | 0.00254 |
| bench_xlarge | 59,049 | 3 | 15.24 | 2.69 | instant | 0.232 | 0.000407 | 0.00251 |
| bench_50x50 | 60,466,176 | 5 | DNF | 30.13 | ≈5 | (~4 min est.) | 0.013 | 1.61 |

## Speedups

**Same algorithm, single-thread (CIBSA full-enum ÷ JuCIB full-enum 1-thread)** —
isolates the pure Julia-vs-Python language gap:

| File | CIBSA / JuCIB-1t |
|---|---:|
| bench_medium | 129× |
| bench_large | 159× |
| bench_typical | 71× |
| bench_xlarge | 66× |

**JuCIB's best vs CIBSA full-enum (exhaustive 10-thread ÷ CIBSA full-enum)** —
algorithm + threading + language combined. These ratios ride on sub-millisecond
JuCIB medians near the timer floor, so treat them as order-of-magnitude
(~10³–10⁴×), not precise:

| File | CIBSA-full / JuCIB-exhaustive |
|---|---:|
| bench_medium | ~1,150× |
| bench_large | ~16,500× |
| bench_typical | ~37,000× |
| bench_xlarge | ~37,000× |

**On the 60 M file:** JuCIB exhaustive (0.013 s, 5 fixed points) is **~400×**
faster than ScenarioWizard (~5 s, same 5) and **~2,300×** faster than CIBSA-MC
(30.1 s, 3 fixed points) — and unlike CIBSA it is exhaustive, so it finds all
five.

## Why JuCIB is fast — the algorithmic speedups

The gap between JuCIB and the other two tools is not just "Julia vs Python."
Three distinct optimisations stack on top of the language difference. Timing the
exhaustive **sweep** at 1 and 10 threads (find-only medians, `test/algo_split.jl`)
separates them on the two 59,049-scenario files:

| | bench_typical | bench_xlarge |
|---|---:|---:|
| Succession-walk, 1 thread (same algorithm as CIBSA) | 0.233 s | 0.220 s |
| Exhaustive fixed-point sweep, 1 thread | 0.00126 s | 0.00121 s |
| Exhaustive fixed-point sweep, 10 threads | 0.000338 s | 0.000278 s |
| **Algorithm alone** (walk → sweep, 1t vs 1t) | **185×** | **182×** |
| **Threading alone** (sweep 1t → 10t) | 3.7× | 4.4× |
| **Combined** | **~690×** | **~790×** |

The pure language gap (CIBSA's Python walk vs JuCIB's Julia walk — same algorithm,
both single-thread) is the 66–159× in the section above; the factors below are
*additional* to that. The headline is that **most of JuCIB's advantage is
algorithmic, not threading**. A fourth optimisation — branch-and-bound pruning —
extends the same idea to spaces too large to enumerate at all (below).

**1. A cheaper question (the big one, ~180×).** The "same algorithm as CIBSA"
finds the kernel by running a global succession from *every* starting scenario —
follow the path scenario → successor → … until it converges or repeats. Each step
recomputes a full impact balance, and many start points walk through the same
intermediates, so the work is well above O(N). JuCIB's `exhaustive=true` sweep
(`_sweep_fixed_points` / `_sweep_chunk!`) instead asks a strictly **local**
question of each scenario: *is it already its own successor?* For each descriptor
it checks whether the currently-selected variant has the maximum impact score; if
every descriptor's current choice is maximal, the scenario is a fixed point. No
path is followed — the kernel is exactly the set of scenarios passing this
one-step test, so a single linear pass replaces N walks.

**2. Early exit.** The test short-circuits inside `_sweep_chunk!`: as soon as
**one** descriptor has a competing variant that beats the current one (strict `>`,
ties keep the incumbent), the scenario can't be a fixed point and the check bails.
Because almost all scenarios are non-fixed (only 5 of 60 M for `bench_50x50`), most
are rejected after scoring just one or two descriptors, long before a full impact
balance is formed.

**3. Allocation-free incremental enumeration.** Instead of decoding each scenario
index with an `inv_signature` divmod and allocating a fresh vector, the sweep walks
the space like an odometer: decode the first scenario of each thread's chunk once,
then **mixed-radix increment** the scenario vector and its CIM row-index vector
in lock-step. The inner loop allocates nothing and never divides.

**4. Embarrassingly parallel.** The space is split into contiguous per-thread chunks
(`Threads.@threads` / `Threads.@spawn`); each thread fills a thread-local kernel
list, concatenated at the end with no de-duplication (every signature is unique).
That adds the 4–8× on top — sub-linear here because the problems are small and the
CPU is hybrid (6 fast P-cores + 8 slower E-cores). On the 60 M file this threading
takes the sweep from **1.18 s (1 thread) to 0.164 s (10 threads)**.

**5. Branch-and-bound for the huge case.** Above 10⁵ scenarios `exhaustive=true`
switches from the sweep to a branch-and-bound search (`_bnb_fixed_points`): it
assigns descriptors depth-first and prunes any subtree whose suffix score bounds
(`_bnb_bounds`) prove no fixed point can lie inside. On `bench_50x50` this visits
**1,361,846 nodes — 2.25 % of the 60 M space** — and finishes in **0.0114 s vs the
full sweep's 0.164 s (≈14×)**, while returning the identical five fixed points.
This is why the 60 M "exhaustive" cell is 0.013 s rather than a fraction of a
second.

### The basin analysis (`find_basins`)

Basin sizes are a different problem: to know *which* attractor each scenario falls
into and tally the counts, you genuinely must follow successions — the cheap local
test and B&B pruning don't apply. v0.2 computes them in **two phases** and sweeps
all 60 M scenarios in **~1.6 s** (down from ~50 s on the v0.1 single-cache walk):

- **Phase 1 — threaded successor table** (`_successor_table!` /
  `_successor_chunk!`). Each thread owns a contiguous chunk and records, for every
  scenario, the signature of its one succession step. The hot operation is the
  impact balance — the element-wise sum of the `ndesc` CIM rows the current
  scenario selects — computed *incrementally*: a mixed-radix odometer advances the
  scenario and the impact-balance buffer is updated by the **row delta** of the one
  digit that changed, so each step touches one CIM column instead of re-summing all
  `ndesc`. The accumulation runs over the transposed CIM (`cim_t`, narrowed to the
  smallest integer type that holds the scores via `_score_type`) so the reads are
  unit-stride, and the loop is `@inbounds @simd`:

  ```julia
  @inbounds for i in 1:ndesc          # seed ib with the selected rows
      r = rows[i]
      @simd for j in 1:ndim
          ib[j] += cimT[j, r]          # contiguous column read → packed AVX adds
      end
  end
  # … then, per odometer step, only the changed digit's row delta:
  @simd for j in 1:ndim
      ib[j] += cimT[j, rnew] - cimT[j, rold]
  end
  ```

  `@inbounds` removes the bounds checks that would otherwise stop the loop
  vectorising; `@simd` lets the compiler emit packed AVX adds; and `ib`, the digit
  vector and row-index vector are pre-allocated once per chunk and reused, so the
  inner loop does **zero heap allocation**. The chunks write disjoint ranges of the
  table, so no synchronisation is needed.

- **Phase 2 — O(N) path-compressed resolve** (`_resolve_and_tally`). A flat `res`
  array records, per scenario, the fixed point (or cycle) it resolves to. Walking
  the successor table from an unresolved scenario marks the chain in-progress; when
  it reaches an already-resolved scenario, a fixed point, or its own chain (a
  cycle), the **whole chain just visited is back-filled** with that result and the
  walk stops. Every scenario is resolved **exactly once** — total work is O(N)
  succession-table hops, not O(N × path length). The final tally over `res` is
  threaded again.

### Net effect

This is why JuCIB's best mode is ~0.4 ms on the mid files where CIBSA's full
enumeration takes 15–16 s, and why the 60 M file goes from "a pure-Python walk that
never finishes" to a complete, verified kernel in **~0.013 s** (branch-and-bound)
plus full basin sizes in **~1.6 s**.

## Verification — do the kernels agree?

**Yes, wherever a full search was run.** Every tool's consistent-scenario set was
reduced to sorted signatures in JuCIB's space and diffed against the
JuCIB-exhaustive reference; every reported scenario was independently confirmed to
satisfy `succession_step(u) == u`. The v0.2 engine returns bit-identical kernels
and basin sizes to v0.1 (both `find_consistent` algorithms and `find_basins` are
cross-checked against each other and a naive oracle in the property tests).

| File | JuCIB (exh) | CIBSA (full) | ScenarioWizard | CIBSA (MC) |
|---|---|---|---|---|
| CIB_global | {13,16,20,21} | = | = | — |
| bench_medium | {191,480} | = | = | — |
| bench_large | {1973,2548} | = | = | — |
| bench_typical | {13785,13839} | = | = | = |
| bench_xlarge | {15921,18086,36255} | = | = | = |
| bench_50x50 | {1984298, 8693546, 11540522, 19807290, 31030843} | DNF | = (all 5) | **subset (3/5)** |

The single discrepancy is **CIBSA in MC mode on the 60 M file**: it recovered only
3 of the 5 fixed points. This is the *expected* behaviour of Monte-Carlo sampling
on a space with tiny basins, **not** a bug — every scenario it *did* report is a
genuine fixed point (no false positives). JuCIB's basin analysis shows why: of the
60,466,176 scenarios, only **119,729 (0.2 %)** converge to any fixed point at all
(the other 60.3 M fall into cycles), and the five basins are very uneven:

| Fixed point (sig) | Basin size | Found by CIBSA-MC? |
|---|---:|:--:|
| 1984298 | 3 | no |
| 19807290 | 95 | yes |
| 8693546 | 578 | yes |
| 31030843 | 5,633 | no |
| 11540522 | 113,420 | yes |

A 10,000-point sample of 60.5 M lands in these basins essentially at random
(expected hits range from 0.0005 to ~19), so *which* fixed points MC finds is
seed-dependent and incomplete. ScenarioWizard and JuCIB-exhaustive both enumerate
the whole space and therefore find all five every time. This is the headline
argument for exhaustive search on large, low-basin problems.

(Tie-breaking was identical across JuCIB and CIBSA — strict `>`, favour current
then lowest index — so no tie-related divergence appeared; SW agreed on every
scenario, so its convention matches too.)

## Caveats

- **JuCIB numbers are v0.2 (2026-07-21); CIBSA/SW are 2026-06-04.** Only JuCIB was
  re-timed. CIBSA and ScenarioWizard are unchanged tools measured on the same
  machine a few weeks earlier, so the cross-tool comparison holds, but the two
  legs were not captured in the same session.
- **CIBSA required a Python-2→3 port to run at all.** The published
  `CIB_sim_anneal.py` is Python-2 source: `print` statements (SyntaxError under
  Py3, so it won't even import), `self.cim[r] = map(...)` (Py3 returns an iterator,
  not a list), a bare `find_consistent()` in the constructor (`NameError`), and a
  mis-named `@kernel.setter` (copy-paste bug → `self.kernel = …` raises). Ten lines
  were changed; the algorithm — succession walk, `find_consistent`, signatures — is
  untouched, so the timings reflect the real CIBSA algorithm. But note the
  published package does not run out of the box on a modern Python.
- **ScenarioWizard timing is coarse.** It is GUI-only; the numbers are
  screenshot-paced wall-clock for the whole operation (parse + search + render),
  not isolated search time. Sub-second results are reported as "instant" because
  GUI timing cannot resolve them. The 60 M number (~5 s) is a 1-call observation,
  not a median.
- **Sub-millisecond JuCIB medians are noisy.** The mid-file `exhaustive` and
  `find_basins` numbers sit near the `time_ns` floor; run-to-run variance of ±50 %
  on a 0.0001–0.0005 s measurement is normal, which is why the "×CIBSA" ratios on
  those rows are quoted as order-of-magnitude.
- **10-thread vs 1-thread.** JuCIB's "exhaustive" and "find_basins" columns use 10
  threads; the only strictly like-for-like single-thread comparison with CIBSA is
  the JuCIB "full-enum 1-thread" column.
- **MC is stochastic.** CIBSA-MC used `seed=999`; the specific fixed points missed
  on the 60 M file would differ with another seed. The *fact* of incompleteness on
  this low-basin problem is robust.
- The "~4 M Cascade Institute" file in the original brief was not available on this
  machine and is out of scope.

## Artifacts & reproduction

Raw results (machine-readable):
- `test/bench_results_julia.json` — JuCIB v0.2, all modes, all files (incl. basin sizes), 2026-07-21.
- `test/bench_results_cibsa.json` — CIBSA full-enum + MC (2026-06-04).
- `test/bench_results_sw.json` — ScenarioWizard kernels (from `.sl`) + observed times (2026-06-04).
- `test/sample_files/<name>_sw.sl` — ScenarioWizard's exported solution sets.

Re-run:
```bash
julia -t 10 --project=. test/bench_threeway.jl   # JuCIB (writes bench_results_julia.json)
julia -t 10 --project=. test/bench_optim.jl      # sweep vs branch-and-bound, B&B node counts
julia -t 1  --project=. test/algo_split.jl       # 1-thread leg of the algorithm/threading split
julia -t 10 --project=. test/algo_split.jl       # 10-thread leg
python test/bench_cibsa.py                        # CIBSA (imports the Py3-ported clone)
julia --project=. test/verify_sw.jl               # verify SW .sl files vs JuCIB
```
ScenarioWizard is GUI-driven: *File ▸ Open* a `.scw`, *Analyse ▸ Consistent
scenarios*, then *Save* the solution set to `<name>_sw.sl`.

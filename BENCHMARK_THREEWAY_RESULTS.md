# Three-way CIB benchmark: ScenarioWizard vs CIBSA vs JuCIB

Wall-clock and correctness comparison of three Cross-Impact Balance (CIB)
consistency-search implementations on a common set of `.scw` files.
Run **2026-06-04**.

- **ScenarioWizard (SW)** — Weimer-Jehle's reference desktop app (GUI, driven via UI automation).
- **CIBSA** — published Python package [`sei-international/cibsa`](https://github.com/sei-international/cibsa) (`CIB_sim_anneal.py`).
- **JuCIB** — `CrossImpactBalances.jl` (this repo).

## TL;DR

1. **All three tools agree on the kernel** for every file where a full search was
   run — identical consistent-scenario sets, down to the signature. The only
   discrepancy is CIBSA in Monte-Carlo mode missing small-basin fixed points on
   the 60 M-scenario file (expected; see Verification).
2. **ScenarioWizard's "much faster than CIBSA" reputation holds** — emphatically.
   SW returns *instantly* on files where CIBSA's pure-Python full enumeration
   takes 15–16 s, and handles the 60 M file in ~5 s where CIBSA cannot enumerate
   it at all.
3. **But JuCIB matches or beats ScenarioWizard.** JuCIB's threaded exhaustive
   search solves the 60 M file in **1.06 s and finds all 5 fixed points** (vs SW's
   ~5 s for the same 5), and is sub-millisecond on the mid-size files. Same
   algorithm, single-thread, JuCIB is **30–95× faster than CIBSA** (pure
   Julia-vs-Python); with its optimised threaded path it is **~2,000–38,000×**
   faster than CIBSA full-enum on the mid files.
4. JuCIB additionally computes full **basin-of-attraction** sizes (e.g. the entire
   60 M-scenario space in ~50 s), which neither of the other tools exposes directly.

## Environment

| Item | Value |
|---|---|
| Machine | Lenovo laptop, Windows 11 Home 10.0.26200 |
| CPU | Intel Core i7-13700H — 14 physical cores (6 P + 8 E) / 20 logical |
| RAM | 15.7 GB |
| Julia | 1.12.5, run with `-t 10` (`Threads.nthreads() = 10`) |
| Python | 3.13.13, numpy 2.4.1, scipy 1.17.0 |
| ScenarioWizard | 5 (file version 1.0.9522.17909), `C:\Program Files (x86)\ScenarioWizard\` |

> Note on threads: the CPU is hybrid (6 performance + 8 efficiency cores). Julia
> ran with 10 threads, the value configured in `JULIA_NUM_THREADS`. Only JuCIB's
> `exhaustive=true` path is multi-threaded.

## Method — what was timed

- **Find-only, median of 3** (warm-up run discarded; `time_ns`/`perf_counter`).
  The model is parsed *once* up front and excluded from the timed region, so JuCIB
  and CIBSA are compared on the search alone. Parse cost is negligible everywhere
  (≤ 0.6 ms even for the 60 M file), so full-load ≈ find-only.
- **JuCIB, three modes:** (a) `find_consistent` single-thread full enumeration
  (`exhaustive=false`, `mc_threshold` above the space size) — the *same algorithm*
  as CIBSA; (b) `find_consistent(…; exhaustive=true)` — the optimised threaded
  fixed-point sweep (JuCIB's best); (c) `find_basins` — full basin analysis.
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
verified as true fixed points.

### CIB_global — 36 scenarios (3 descriptors)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.000017 | 4 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000077 | 4 | ✓ |
| JuCIB | find_basins | 1 | 0.000004 | 4 | ✓ |
| CIBSA | full-enum | 1 | 0.001626 | 4 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 4 | ✓ |

Kernel = `{13, 16, 20, 21}`.

### bench_medium — 1,024 scenarios (5×4)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.00196 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000076 | 2 | ✓ |
| JuCIB | find_basins | 1 | 0.000205 | 2 | ✓ |
| CIBSA | full-enum | 1 | 0.151 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{191, 480}`.

### bench_large — 6,561 scenarios (8×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.0282 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000070 | 2 | ✓ |
| JuCIB | find_basins | 1 | 0.00196 | 2 | ✓ |
| CIBSA | full-enum | 1 | 2.668 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{1973, 2548}`.

### bench_typical — 59,049 scenarios (10×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.476 | 2 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000754 | 2 | ✓ |
| JuCIB | find_basins | 1 | 0.0285 | 2 | ✓ |
| CIBSA | full-enum | 1 | 16.35 | 2 | ✓ |
| CIBSA | MC (10 k) | 1 | 2.729 | 2 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 2 | ✓ |

Kernel = `{13785, 13839}`.

### bench_xlarge — 59,049 scenarios (10×3)

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | full-enum | 1 | 0.474 | 3 | ✓ |
| JuCIB | **exhaustive** | 10 | 0.000859 | 3 | ✓ |
| JuCIB | find_basins | 1 | 0.0279 | 3 | ✓ |
| CIBSA | full-enum | 1 | 15.24 | 3 | ✓ |
| CIBSA | MC (10 k) | 1 | 2.686 | 3 | ✓ |
| ScenarioWizard | GUI full | — | instant¹ | 3 | ✓ |

Kernel = `{15921, 18086, 36255}`.

### bench_50x50 — 60,466,176 scenarios (15 descriptors) — the stress test

| Tool | Mode | Threads | Time (s) | Kernel | Match |
|---|---|---:|---:|---:|:--:|
| JuCIB | **exhaustive** | 10 | **1.060** | 5 | ✓ (reference) |
| JuCIB | find_basins | 1 | 50.43 | 5 | ✓ |
| JuCIB | full-enum | 1 | not run² | — | — |
| CIBSA | full-enum | 1 | **DNF**³ | — | — |
| CIBSA | MC (10 k) | 1 | 30.13 | **3** | partial⁴ |
| ScenarioWizard | GUI full | — | **≈ 2–6 s**⁵ | 5 | ✓ |

Kernel (exhaustive) = `{1984298, 8693546, 11540522, 19807290, 31030843}`.

¹ Result appeared instantly with no progress dialog → below the ~0.5 s resolution
of GUI timing.
² Single-thread succession-walk enumeration of 60 M scenarios was skipped; naive
extrapolation from `bench_typical` is on the order of ~8 min. JuCIB's threaded
exhaustive path (1.06 s) is the real number.
³ Pure-Python full enumeration of 60 M scenarios is infeasible (hours) — not run,
per the benchmark's "abandon after 30 min" rule.
⁴ CIBSA-MC (seed 999, 10 k samples) found 3 of the 5 fixed points. See Verification.
⁵ Coarse — a "Computing…" progress dialog appeared and the result rendered within
the next ~5 s.

## Cross-file summary (find-only, seconds)

| File | Scenarios | Kernel | CIBSA full | CIBSA MC | SW | JuCIB 1-thread | JuCIB exhaustive (10t) | JuCIB basins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CIB_global | 36 | 4 | 0.0016 | — | instant | 0.000017 | 0.000077 | 0.000004 |
| bench_medium | 1,024 | 2 | 0.151 | — | instant | 0.00196 | 0.000076 | 0.000205 |
| bench_large | 6,561 | 2 | 2.668 | — | instant | 0.0282 | 0.000070 | 0.00196 |
| bench_typical | 59,049 | 2 | 16.35 | 2.73 | instant | 0.476 | 0.000754 | 0.0285 |
| bench_xlarge | 59,049 | 3 | 15.24 | 2.69 | instant | 0.474 | 0.000859 | 0.0279 |
| bench_50x50 | 60,466,176 | 5 | DNF | 30.13 | ≈5 | (~8 min est.) | 1.060 | 50.43 |

## Speedups

**Same algorithm, single-thread (CIBSA full-enum ÷ JuCIB full-enum 1-thread)** —
isolates the pure Julia-vs-Python language gap:

| File | CIBSA / JuCIB-1t |
|---|---:|
| bench_medium | 77× |
| bench_large | 95× |
| bench_typical | 34× |
| bench_xlarge | 32× |

**JuCIB's best vs CIBSA full-enum (exhaustive 10-thread ÷ CIBSA full-enum)** —
algorithm + threading + language combined:

| File | CIBSA-full / JuCIB-exhaustive |
|---|---:|
| bench_medium | ~1,990× |
| bench_large | ~38,000× |
| bench_typical | ~21,700× |
| bench_xlarge | ~17,700× |

**On the 60 M file:** JuCIB exhaustive (1.06 s, 5 fixed points) is ~5× faster than
ScenarioWizard (~5 s, same 5) and ~28× faster than CIBSA-MC (30.1 s, 3 fixed
points) — and unlike CIBSA it is exhaustive, so it finds all five.

## Why JuCIB is fast — the algorithmic speedups

The gap between JuCIB and the other two tools is not just "Julia vs Python."
Three distinct optimisations stack on top of the language difference. Timing the
exhaustive sweep at 1 and 10 threads (find-only medians, `test/algo_split.jl`)
separates them on the two 59,049-scenario files:

| | bench_typical | bench_xlarge |
|---|---:|---:|
| Succession-walk, 1 thread (same algorithm as CIBSA) | 0.51 s | 0.41 s |
| Exhaustive fixed-point sweep, 1 thread | 0.0019 s | 0.0027 s |
| Exhaustive fixed-point sweep, 10 threads | 0.0005 s | 0.0003 s |
| **Algorithm alone** (walk → sweep, 1t vs 1t) | **267×** | **148×** |
| **Threading alone** (sweep 1t → 10t) | 3.8× | 8.0× |
| **Combined** | **~1020×** | **~1190×** |

The pure language gap (CIBSA's Python walk vs JuCIB's Julia walk — same algorithm,
both single-thread) is the 32–95× in the section above; the factors below are
*additional* to that. The headline is that **most of JuCIB's advantage is
algorithmic, not threading**.

**1. A cheaper question (the big one, ~150–270×).** The "same algorithm as CIBSA"
finds the kernel by running a global succession from *every* starting scenario —
follow the path scenario → successor → … until it converges or repeats. Each step
recomputes a full impact balance, and many start points walk through the same
intermediates, so the work is well above O(N). JuCIB's `exhaustive=true` path
instead asks a strictly **local** question of each scenario: *is it already its own
successor?* For each descriptor it checks whether the currently-selected variant
has the maximum impact score; if every descriptor's current choice is maximal, the
scenario is a fixed point. No path is followed — the kernel is exactly the set of
scenarios passing this one-step test, so a single linear pass replaces N walks.

**2. Early exit.** The test short-circuits (`src/CrossImpactBalances.jl` L401–428):
as soon as **one** descriptor has a competing variant that beats the current one,
the scenario can't be a fixed point and the check bails. Because almost all
scenarios are non-fixed (only 5 of 60 M for `bench_50x50`), most are rejected after
scoring just one or two descriptors, long before a full impact balance is formed.

**3. Allocation-free incremental enumeration.** Instead of decoding each scenario
index with an `inv_signature` divmod and allocating a fresh vector, the sweep walks
the space like an odometer: decode the first scenario of each thread's chunk once,
then **mixed-radix increment** the scenario vector `v` and its CIM row-index vector
`tndx` in lock-step (L434–445). The inner loop allocates nothing and never divides.

**4. Embarrassingly parallel.** The space is split into contiguous per-thread chunks
(`Threads.@threads`); each thread fills a thread-local kernel list, concatenated at
the end with no de-duplication (every signature is unique). That adds the 4–8× on
top — sub-linear here because the problems are small and the CPU is hybrid (6 fast
P-cores + 8 slower E-cores). On the 60 M file this threading takes the sweep from
3.67 s (1 thread) to ~0.6–1.1 s (10 threads).

### SIMD: the basin analysis (`find_basins`)

Basin sizes are a different problem: to know *which* attractor each scenario falls
into and tally the counts, you genuinely must follow successions — the cheap local
test doesn't apply. Two things keep `find_basins` fast enough to sweep all 60 M
scenarios in ~50 s single-threaded:

- **O(N) memoised succession.** A flat `cache` array records, per scenario, the
  fixed point it resolves to (or "cycle"). When a chain reaches an already-resolved
  scenario, the whole chain just visited is back-filled with that result and the
  walk stops (L537–544). Every scenario is resolved **exactly once** — total work is
  O(N) succession steps, not O(N × path length).

- **A vectorised, allocation-free inner kernel.** The hot operation is the impact
  balance: the element-wise sum of the `ndesc` CIM rows the current scenario
  selects. Rather than call the generic `impact_balance` (which allocates a fresh
  vector and goes through index indirection), `find_basins` inlines it as a
  row-at-a-time accumulation into a reused buffer, annotated `@inbounds @simd`
  (L506–514):

  ```julia
  r1 = tndx[1]
  @inbounds @simd for j in 1:ndim          # seed ib with the first selected row
      ib[j] = cim[r1, j]
  end
  for ki in 2:ndesc                         # add each remaining selected row
      r = @inbounds tndx[ki]
      @inbounds @simd for j in 1:ndim
          ib[j] += cim[r, j]
      end
  end
  ```

  `@inbounds` removes the bounds checks that would otherwise stop the loop
  vectorising; `@simd` lets the compiler emit packed AVX adds, so several
  `ib[j] += …` lanes are computed per instruction; and `ib`, `w`, `tndx`, `history`
  are pre-allocated once and reused, so the inner loop does **zero heap allocation**.
  Picking the best variant per descriptor (L516–527) runs in the same `@inbounds`
  block. (Since each accumulation runs along a CIM *row* while Julia stores matrices
  column-major, a further refinement is to keep the CIM transposed so those reads
  are unit-stride — the direction the `cim-rowmajor` branch name points at.)

### Net effect

This is why JuCIB's best mode is ~0.3–0.5 ms on the mid files where CIBSA's full
enumeration takes 15–16 s, and why the 60 M file goes from "a pure-Python walk that
never finishes" to a complete, verified kernel in ~0.6–1.1 s plus full basin sizes
in ~50 s.

## Verification — do the kernels agree?

**Yes, wherever a full search was run.** Every tool's consistent-scenario set was
reduced to sorted signatures in JuCIB's space and diffed against the
JuCIB-exhaustive reference; every reported scenario was independently confirmed to
satisfy `succession_step(u) == u`.

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

- **CIBSA required a Python-2→3 port to run at all.** The published
  `CIB_sim_anneal.py` is Python-2 source: `print` statements (SyntaxError under
  Py3, so it won't even import), `self.cim[r] = map(...)` (Py3 returns an iterator,
  not a list), a bare `find_consistent()` in the constructor (`NameError`), and a
  mis-named `@kernel.setter` (copy-paste bug → `self.kernel = …` raises). Ten lines
  were changed (`git diff` in `D:\GitHub\cibsa`); the algorithm — succession walk,
  `find_consistent`, signatures — is untouched, so the timings reflect the real
  CIBSA algorithm. But note the published package does not run out of the box on a
  modern Python.
- **ScenarioWizard timing is coarse.** It is GUI-only; the numbers are
  screenshot-paced wall-clock for the whole operation (parse + search + render),
  not isolated search time. Sub-second results are reported as "instant" because
  GUI timing cannot resolve them. The 60 M number (~5 s) is a 1-call observation,
  not a median.
- **10-thread vs 1-thread.** JuCIB's "exhaustive" column uses 10 threads; the only
  strictly like-for-like single-thread comparison with CIBSA is the JuCIB
  "full-enum 1-thread" column.
- **MC is stochastic.** CIBSA-MC used `seed=999`; the specific fixed points missed
  on the 60 M file would differ with another seed. The *fact* of incompleteness on
  this low-basin problem is robust.
- The "~4 M Cascade Institute" file in the original brief was not available on this
  machine and is out of scope.

## Artifacts & reproduction

Raw results (machine-readable):
- `test/bench_results_julia.json` — JuCIB, all modes, all files (incl. basin sizes).
- `test/bench_results_cibsa.json` — CIBSA full-enum + MC.
- `test/bench_results_sw.json` — ScenarioWizard kernels (from `.sl`) + observed times.
- `test/sample_files/<name>_sw.sl` — ScenarioWizard's exported solution sets.

Re-run:
```bash
julia -t 10 --project=. test/bench_threeway.jl   # JuCIB
python test/bench_cibsa.py                        # CIBSA (imports the Py3-ported clone)
julia --project=. test/verify_sw.jl               # verify SW .sl files vs JuCIB
```
ScenarioWizard is GUI-driven: *File ▸ Open* a `.scw`, *Analyse ▸ Consistent
scenarios*, then *Save* the solution set to `<name>_sw.sl`.

---
title: 'CrossImpactBalances.jl: fast, exact cross-impact balance scenario analysis in Julia'
tags:
  - Julia
  - Python
  - scenario analysis
  - cross-impact balance
  - foresight
  - qualitative systems analysis
authors:
  - name: Iain Morrow
    # TODO(author): replace with your real ORCID before submission — JOSS requires it.
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
affiliations:
  # TODO(author): confirm affiliation wording (institution, or "Independent Researcher, <country>").
  - name: Independent Researcher, United Kingdom
    index: 1
date: 29 July 2026
bibliography: paper.bib
---

# Summary

Cross-impact balance (CIB) analysis [@weimerjehle2006] is an established method
for constructing qualitative scenarios of systems that resist quantitative
modelling, such as energy transitions, socio-technical change, and public
health. Experts describe a system as a set of descriptors, each taking one of a
few discrete states, and judge how each state of each descriptor promotes or
inhibits each state of every other. A scenario — one state per descriptor — is
*internally consistent* if no descriptor's chosen state is beaten, under the
combined influences of all the others, by an alternative state. The set of
consistent scenarios (the *kernel*) comprises the futures that the judgement
network can support.

`CrossImpactBalances.jl` is a Julia package, with no dependencies outside the
standard library, that computes CIB kernels *exactly* — no sampling and no
completeness caveat — on scenario spaces of up to hundreds of millions of
combinations, in milliseconds to seconds on a laptop. It also computes complete
basin-of-attraction statistics (which attractor every scenario converges to,
and how large each basin is), an analysis that no existing CIB tool offers. It
reads the `.scw` files of the reference tool ScenarioWizard, and is usable from
Julia, from Python (in-process, or through a compiled shared library that
requires no Julia installation), and — through a local-browser desktop
application with a Windows installer — by non-programmers.

# Statement of need

CIB underpins a growing body of scenario research, including
internal-consistency screening of the IPCC SRES storylines [@schweizer2012] and
systematic construction of shared socioeconomic pathway elements
[@schweizer2014]; the method and its applications are surveyed by
@weimerjehle2023. A CIB matrix with $d$ descriptors of $v$ states each defines
$v^d$ scenarios, so realistic matrices quickly outgrow what exhaustive analysis
has conventionally been thought to permit: fifteen descriptors of three to four
states already yield tens of millions of scenarios. Three needs of working
scenario researchers are unmet by existing tools.

**Exactness at scale.** Above modest problem sizes, existing programmable
tooling falls back on Monte-Carlo sampling or simulated annealing, and sampling
provably misses consistent scenarios whose basins of attraction are small. This
is not hypothetical: on a benchmark 60-million-scenario matrix, only 0.2% of
scenarios converge to any fixed point, basin sizes span 3 to 113,420, and a
10,000-point Monte-Carlo run recovers only three of the five consistent
scenarios. Because the kernel *is* the result of a CIB study, a silently
incomplete kernel is a wrong answer. `CrossImpactBalances.jl` is always
exhaustive.

**Throughput for robustness and uncertainty analysis.** Cross-impact
judgements are elicited, ordinal, and contestable, so good practice varies
them: sensitivity sweeps over individual judgements, batch runs over ensembles
of matrices, and uncertainty analysis over judgement perturbations all require
thousands of solves. Millisecond kernel searches, in-place matrix editing
without re-parsing, and batch/sweep helpers in the Python wrapper make such
meta-analyses interactive rather than overnight jobs.

**Dynamics beyond the kernel.** CIB defines a deterministic succession
operator whose fixed points are the consistent scenarios. The basin of
attraction of each fixed point — how much of the scenario space converges to it
— is a natural robustness measure, and cycles of the operator are of
methodological interest in their own right, yet no CIB tool exposes either.
`find_basins` returns exact basin sizes and cycle counts for the entire space.
The succession rule itself is also pluggable, supporting research on CIB
dynamics: a user-defined rule needs one method, and rules that declare a
fixed-point margin inherit the fast search paths.

# State of the field

Two implementations dominate CIB practice. ScenarioWizard [@scenariowizard],
the method author's free reference implementation, is a closed-source Windows
GUI with no programmatic API and no batch mode. `cibsa` [@cibsa], an
open-source Python module from the Stockholm Environment Institute, is the
direct ancestor of this package; it is exact only on small spaces and samples
(Monte-Carlo or simulated annealing) on large ones. Table 1 summarises a
three-way comparison on a common set of `.scw` files (Intel i7-13700H; Julia
with 10 threads; ScenarioWizard timings are coarse wall-clock observations of
the GUI; reproduction scripts ship in the repository).

| Problem | Scenarios | `cibsa` full | `cibsa` MC (10k) | ScenarioWizard | `CrossImpactBalances.jl` |
|---|---:|---:|---:|---:|---:|
| `bench_typical` | 59,049 | 16.4 s | 2.7 s | < 1 s | 0.4 ms |
| `bench_50x50` | 60,466,176 | infeasible | 30.1 s, finds 3/5 | ~5 s, finds 5/5 | **13 ms, finds 5/5** |

: Consistency-search wall-clock times (find-only medians). The headline gap on
`bench_typical` — 16.4 s versus 0.4 ms, roughly 37,000× — combines the replaced
search algorithm with 10-way threading; like-for-like (the same
succession-walk algorithm, both single-threaded), the Julia engine is 66–159×
faster than `cibsa` across the benchmark files. On the 60-million-scenario
stress file the package also computes the full basin-of-attraction analysis in
1.6 s — an analysis the other tools do not offer at any size.

Decomposing the gap in Table 1 shows that most of it is methodological rather
than incidental. Two factors are ordinary engineering: expressing the same
succession-walk search in Julia [@bezanson2017] rather than pure Python
accounts for the 66–159×, and multi-threading adds only a further ~4× on the
ten-thread benchmark machine. The remaining orders of magnitude come from
three algorithmic changes, none present in `cibsa` or, to our knowledge, any
other CIB tool:

- **A local consistency test instead of succession walks** (~180× at equal
  thread counts). `cibsa` finds the kernel by running a succession walk from
  every scenario, each step recomputing a full impact balance and many walks
  re-visiting the same intermediate states. This package instead asks each
  scenario a strictly local question — is any descriptor's chosen state beaten
  by an alternative? — and abandons the scenario at the first beaten
  descriptor, so a single linear pass replaces $v^d$ walks.
- **Incremental SIMD scoring.** Scenarios are enumerated by a mixed-radix
  odometer so that consecutive scenarios differ in one descriptor state, and
  the impact balance is updated by that state's score-row delta — a vectorised
  loop over narrow integers with no allocation — instead of being recomputed
  from scratch each time.
- **Exact branch-and-bound** (a further ~50× over the sweep on the
  60-million-scenario file, with fundamentally better scaling). Descriptors
  are assigned depth-first, and precomputed suffix bounds on each rival
  variant's score *difference* — taken over the undecided descriptors'
  shared variant choices, which is strictly tighter than bounding the two
  scores independently — prune any subtree that provably contains no
  consistent scenario; on the stress file the search visits 0.6% of the
  space. Because its cost scales with the
  pruned tree rather than the scenario space, it extends exact analysis
  beyond what enumeration can reach; a visited-node budget falls back to the
  sweep on weakly-coupled matrices where pruning cannot pay off, preserving
  exactness.

Basin-of-attraction analysis is not a speedup but a new capability, with its
own two-phase design (see Software design). These changes replace `cibsa`'s
core rather than tune it — which is why the work could not have been
contributed incrementally upstream — and the package also pursues goals
outside `cibsa`'s scope: an embeddable C ABI, source-free Python
distribution, and desktop packaging. The Python wrapper nevertheless
preserves a `cibsa`-adjacent workflow so existing users can migrate with
little friction.

Correctness is established two ways: a three-way validation showing identical
kernels — down to individual scenario signatures, with every reported scenario
independently re-verified as a true fixed point — across ScenarioWizard,
`cibsa`, and this package on the benchmark corpus; and roughly 1,340
property-based test assertions that pin every engine against brute-force
oracles on randomised instances.

# Software design

The engine is pure Julia with zero package dependencies, which keeps it
auditable, trivially installable, and compilable into a standalone shared
library. The sweep (chosen automatically below $10^5$ scenarios), the
branch-and-bound search (chosen above it), and the basin engine return
bit-identical results. The basin analysis runs in two phases: a threaded pass
fills a flat successor table, reusing the incremental SIMD scoring described
above, and a path-compressed resolution pass then assigns every scenario to
its fixed point or cycle in O(N) total work with O(1) cycle detection. Memory
is independent of thread count (~8 bytes per scenario), which brought a
408-million-scenario model from an infeasible ~13 GB under a per-thread
design to 1.6 GB.

Succession dynamics is an extension point: a new rule subtypes
`SuccessionRule` and defines a single `succession_step` method, and every
analysis routine accepts it; declaring a fixed-point margin lets a rule inherit
the sweep and branch-and-bound searches. A Gauss–Seidel-style
`SequentialSuccession` ships as a second built-in rule.

Around the engine sit a Python package (running the engine in-process via
`juliacall`, or driving a PackageCompiler-built `libcib` through `ctypes` with
no Julia installation), a JSON-based C ABI for embedding from other languages,
a local-browser desktop application packaged as a Windows MSI installer for
non-programming practitioners, and a documentation site with worked examples.

# Research impact statement

The package makes exact kernels and complete basin statistics routine at
scales ($10^7$–$10^8$ scenarios) where current practice either samples — with
the demonstrated risk of silently missing consistent scenarios — or is
confined to a GUI without programmatic access. The benchmark corpus and
scripts, which to our knowledge constitute the first published cross-validation
of CIB implementations (ScenarioWizard, `cibsa`, and this package agreeing on
every kernel), are included in the repository and are reusable for validating
future tools.
<!-- TODO(author): add 1–3 sentences on your own concrete research use of the
package (project, matrix sizes, findings, any collaborators/users). JOSS
expects demonstrated use at minimum by the developers; do not leave this
placeholder in the submitted version. -->
The desktop application and source-free distribution target the practitioner
community that currently depends on ScenarioWizard, lowering the barrier to
scripted, reproducible CIB studies.

# AI usage disclosure

Generative AI tools (Anthropic's Claude) were used during the development of
the software and in the drafting of this manuscript, both as a coding
assistant and as a source of algorithmic suggestions — notably the
branch-and-bound search strategy, which was proposed by the AI and then
reviewed, adapted, and validated by the author. All algorithms, benchmarks,
and text were reviewed and validated by the author; correctness of the
software is established independently of any AI involvement by property-based
tests against brute-force oracles and by cross-validation against
ScenarioWizard and `cibsa`. The author takes full responsibility for the
software and this paper.

# Acknowledgements

This package began as a reimplementation of the Stockholm Environment
Institute's open-source `cibsa` [@cibsa]; its MIT-licensed code informed the
original engine, and its copyright notice is preserved. Wolfgang Weimer-Jehle's
CIB method and ScenarioWizard software define the `.scw` file format and
provided the reference results used for validation.

# References

# Contributing to CrossImpactBalances.jl

Thanks for your interest in improving the package. Contributions of all kinds
are welcome — bug reports, documentation fixes, new succession rules,
performance work, and example models.

## Reporting problems

Open an issue at
[github.com/IainDM/CrossImpactBalances/issues](https://github.com/IainDM/CrossImpactBalances/issues).
For bugs, please include the Julia (or Python) version, the smallest `.scw`
file or matrix that reproduces the problem, and the exact call and output.
Incorrect results are treated as the highest-priority bugs — the package's
core promise is exactness.

## Seeking support

Usage questions are welcome as GitHub issues too — label them as questions.
The [documentation site](https://IainDM.github.io/CrossImpactBalances.jl/) and
the worked examples in [`examples/`](examples/) are the best starting points.

## Contributing changes

1. Fork the repository and create a branch from `main`.
2. Make your change. For engine changes, please keep results exact: every
   search algorithm is pinned against brute-force oracles in
   `test/property_tests.jl`, and new algorithms should be added to those
   cross-checks.
3. Run the tests: `julia --project=. -e 'using Pkg; Pkg.test()'` (and
   `pytest` from `python/` if you touched the Python wrapper).
4. Open a pull request describing what changed and why. Benchmarks
   (`test/bench_optim.jl`, `test/bench_threeway.jl`) are appreciated for
   performance claims.

If you touch the search engine and you have Weimer-Jehle's corpus staged (see
[`test/WWJ_CORPUS.md`](test/WWJ_CORPUS.md)), also run
`julia -t auto --project=. test/verify_wwj.jl`. It checks 38 real ScenarioWizard
models — including 19 against ScenarioWizard's own answers, at up to 10²⁴
scenarios — and exits non-zero on any disagreement. A subset is part of
`Pkg.test()` when the corpus is present, and skips when it is not, so a green
suite on a fresh clone does not mean it ran. `JUCIB_WWJ_FULL=1` adds the two
slowest paired models and the ones with no reference solution; `Pkg.test()`
runs with `--check-bounds=yes`, which makes those searches roughly 7× dearer
than the standalone verifier does.

Third-party model collections are welcome, but they go in their own directory
with a provenance document recording who supplied them and on what terms, rather
than into `test/sample_files/` alongside our own fixtures.

New succession rules are a particularly easy first contribution: subtype
`SuccessionRule`, define one `succession_step` method, and (optionally) a
`fixed_point_margin` — see
[`examples/04_custom_succession_rule.jl`](examples/04_custom_succession_rule.jl).

## Code of conduct

Be respectful and constructive. Reports of unacceptable behaviour can be sent
privately to the maintainer.

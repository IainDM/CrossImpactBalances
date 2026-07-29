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

New succession rules are a particularly easy first contribution: subtype
`SuccessionRule`, define one `succession_step` method, and (optionally) a
`fixed_point_margin` — see
[`examples/04_custom_succession_rule.jl`](examples/04_custom_succession_rule.jl).

## Code of conduct

Be respectful and constructive. Reports of unacceptable behaviour can be sent
privately to the maintainer.

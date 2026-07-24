# Benchmarks

Performance harness for `find_basins` (exhaustive basin-of-attraction analysis).
All scripts run against the sample models in `../test/sample_files` and take an optional
model argument (`nested`, `50x50`, `typical`). Run with the package project active and the
desired thread count:

```
julia --project=. -t 8 bench/<script>.jl [model]
```

| script | what it measures |
|---|---|
| [`stage_breakdown.jl`](stage_breakdown.jl) | per-stage wall time — parallel successor-table build vs resolve+tally. Sweep `-t 1 2 4 8` to see scaling. |
| [`determinism.jl`](determinism.jl) | runs `find_basins` N times and asserts identical output; prints a fingerprint that must match across thread counts (race/nondeterminism check). |
| [`ceilings.jl`](ceilings.jl) | machine STREAM bandwidth + random-chase latency, for interpreting the numbers as compute- vs memory-bound. |

See also [`../test/benchmark_basins.jl`](../test/benchmark_basins.jl) for the end-to-end
timing + correctness benchmark against pinned references.

## Reference numbers

On an Intel i7-9700 (8c/no-SMT, dual-channel DDR4-2667, AVX2), Julia 1.12, 8 threads:
`CIB_nested` (408,146,688 states) ≈ **3.0 s**; `bench_50x50` (60,466,176) ≈ **0.47 s**.
`find_basins` on these models is compute-bound — the table build (SIMD argmax) dominates,
the resolve is cache-served, and nothing is DRAM-bandwidth or -latency bound.

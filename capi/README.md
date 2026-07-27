# capi — C-callable shared library

This directory compiles the CrossImpactBalances engine into a **standalone
shared library** (`libcib.so` / `.dylib` / `.dll`) so it can be driven from
non-Julia callers — notably the Python package's *native* backend — **without
shipping any Julia source**. The library contains machine code plus a bundled
Julia runtime.

## What's here

| File | Purpose |
|------|---------|
| `src/CIBCApi.jl` | The C ABI: `@ccallable`, JSON-in/JSON-out functions over integer model handles. |
| `src/json.jl` | Vendored dependency-free JSON parser/serializer (copy of `mcp/json.jl`). |
| `Project.toml` | The `CIBCApi` package; depends on `CrossImpactBalances`. |
| `precompile.jl` | Execution trace so PackageCompiler bakes the hot paths into the library. |

## Building

From the repository root, on the platform you want to target:

```bash
julia --project=build build/build_library.jl [OUTDIR]
```

Default `OUTDIR` is `build/cib-lib/`, containing:

```
cib-lib/
  lib/libcib.{so,dylib,dll}   compiled engine (+ Julia runtime .so/.dylib/.dll)
  include/cib.h, julia_init.h C headers
  share/julia/…               runtime support
```

The bundle is relocatable and needs no Julia installation to run. The build is
**OS/architecture-native**: run it on each platform you ship
(Linux/macOS/Windows × x86_64/arm64) to produce that platform's library.

## The C ABI (for reference)

Runtime lifecycle (from `julia_init.h`): call `init_julia(argc, argv)` once
before any `cib_*` call and `shutdown_julia(retcode)` at the end. Pass
`-t auto`/`-tN` in `argv` to enable engine threads.

Every `cib_*` function takes/returns JSON as C strings (`char *`). Returned
strings are heap-allocated and must be released with `cib_free_string`. Models
live behind integer handles; release them with `cib_free`. Scenarios are JSON
arrays of 0-based variant indices.

| Function | JSON in → out |
|----------|---------------|
| `cib_load(path)` | → `{handle, descriptors, variants, n_descriptors, n_scenarios}` |
| `cib_consistent(h, opts)` | `{rule?, algorithm?}` → `{scenarios}` |
| `cib_basins(h, opts)` | `{rule?}` → `{fixed_points, basin_sizes, cycle_count, total}` |
| `cib_impact_balance(h, scenario)` | `[i,…]` → `{ib}` |
| `cib_succession(h, req)` | `{scenario, rule?, max_steps?}` → `{steps, cycle_length}` |
| `cib_signature(h, scenario)` | `[i,…]` → `{signature}` |
| `cib_inv_signature(h, sig)` | `s` → `{scenario}` |
| `cib_set_impact(h, req)` | `{source:[d,v], target:[d,v], value}` → `{old}` |
| `cib_get_impact(h, req)` | `{source:[d,v], target:[d,v]}` → `{value}` |
| `cib_matrix(h)` | → `{matrix}` |
| `cib_copy(h)` | → `{handle}` |
| `cib_free(h)` / `cib_free_string(ptr)` | void |

Errors never cross the boundary as exceptions: a failed call returns
`{"ok": false, "error": "…"}`.

`cib_consistent` previously also accepted `exhaustive`, `ignore_cycles` and
`seed`. Those are gone along with the simulated annealing and Monte-Carlo
sampling they configured: the search now always covers the whole scenario
space and is fully deterministic, so there is nothing to seed and no partial
search to opt out of. Unknown keys are ignored, so old callers keep working.

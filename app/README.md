# Cross-Impact Balances — desktop app

A small local application that wraps
[`CrossImpactBalances.jl`](../) in a point-and-click UI:

1. **Browse** to a ScenarioWizard `.scw` file.
2. **Find Consistent Scenarios** — lists every consistent scenario (fixed
   point of the succession map), fast even for very large models
   (branch-and-bound).
3. **Find Basins** — full basin-of-attraction analysis: each consistent
   scenario with the size of its basin (how many starting scenarios converge
   to it), the share of the space that converges to a fixed point, and the
   number of starts that fall into non-fixed-point cycles.
4. **Export basin CSV** — the full basin table (all descriptor columns, basin
   sizes and fractions) for use in Excel, R, or pandas.

The app runs entirely on your machine. It serves a single page to your
browser over `localhost` — nothing is uploaded anywhere. It has **no external
Julia dependencies** (just the standard library plus the analysis engine), so
it packages cleanly into a self-contained installer.

## Run from source

From the repository root:

```bash
julia --project=app -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=app -t auto app/run.jl
```

This opens `http://127.0.0.1:8071/` in your default browser. Use `-t auto`
so the exhaustive search and basin analysis use all CPU cores. Click **Quit**
in the page (or close the terminal window) to stop the server.

## Build a Windows installer (.msi)

See [`../build/README.md`](../build/README.md). The MSI compiles this app into
a standalone executable with PackageCompiler and wraps it with WiX, so it can
be handed to someone without installing Julia or the source code. A GitHub
Actions workflow builds the `.msi` automatically on a Windows runner.

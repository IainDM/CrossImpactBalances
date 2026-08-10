# Cross-Impact Balances — desktop app

A small local application that wraps
[`CrossImpactBalances.jl`](../) in a point-and-click UI:

1. **Browse** to a ScenarioWizard `.scw` file.
2. **Find Consistent Scenarios** — lists every consistent scenario (fixed
   point of the succession map), fast even for very large models
   (branch-and-bound). A model with more than a couple of thousand of them is
   reported by its exact kernel size instead of listed — the count is reached
   through the influence map, so it arrives in seconds even where the kernel
   itself runs to billions of scenarios and could never be held in memory.
3. **Find Basins** — full basin-of-attraction analysis: each consistent
   scenario with the size of its basin (how many starting scenarios converge
   to it), the share of the space that converges to a fixed point, and the
   number of starts that fall into non-fixed-point cycles.
4. **Export basin CSV** — the full basin table (all descriptor columns, basin
   sizes and fractions) for use in Excel, R, or pandas. A spreadsheet holds far
   more rows than a page, so a basin analysis with too many consistent
   scenarios to draw is still written to CSV — up to a million of them, past
   which only the size is reported.

The app runs entirely on your machine. It serves a single page to your
browser over `localhost` — nothing is uploaded anywhere. It has **no external
Julia dependencies** (just the standard library plus the analysis engine), so
it packages cleanly into a self-contained installer.

## Run from source

From the repository root:

```bash
# first time only — resolve the engine from the repo root and install deps
julia --project=app -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=app -t auto,1 app/run.jl
```

This opens `http://127.0.0.1:8071/` in your default browser. `auto,1` gives the
exhaustive search and basin analysis all CPU cores, plus one interactive thread
— that thread is what keeps the page answering, and **Quit** working, while a
long analysis has every other core busy. Click **Quit** in the page (or close
the terminal window) to stop the server.

## Build a Windows installer (.msi)

See [`../build/README.md`](../build/README.md). The MSI compiles this app into
a standalone executable with PackageCompiler and wraps it with WiX, so it can
be handed to someone without installing Julia or the source code. A GitHub
Actions workflow builds the `.msi` automatically on a Windows runner.

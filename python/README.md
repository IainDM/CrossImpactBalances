# crossimpactbalances (Python)

A Pythonic interface to
[CrossImpactBalances.jl](https://github.com/IainDM/CrossImpactBalances), the
fast Julia engine for Cross-Impact Balance (CIB) scenario analysis. The Julia
engine runs **in-process** via [`juliacall`](https://github.com/JuliaPy/PythonCall.jl),
so a parsed model stays resident: you can run many analyses, batch over many
`.scw` files, and tweak individual cross-impact values in place — all without
re-parsing.

## Install

From the repository root (editable install links against the in-repo Julia
source):

```bash
pip install -e "python/[pandas]"
```

The first call into the engine provisions a private Julia and develops the
local `CrossImpactBalances.jl` package (via [juliapkg](https://github.com/JuliaPy/pyjuliapkg));
this is a one-time cost.

### Multi-threading

The exhaustive-search and basin routines parallelise across Julia threads. Set
this **before** the first engine call:

```bash
export PYTHON_JULIACALL_THREADS=auto
```

## Quick start

```python
from crossimpactbalances import Model, run_models, sweep_impact

# Load a ScenarioWizard .scw model (parsed once, kept resident).
m = Model.load("test/sample_files/CIB_global.scw")

m.descriptors                      # ['WTRD', 'WSEC', 'WECO']
m.variants                         # {'WTRD': ['FT', 'Ntl', 'Mix'], ...}
m.n_scenarios                      # 36

# Consistent scenarios come back as {descriptor: variant} name maps.
for sc in m.consistent_scenarios():
    print(m.signature(sc), sc)

# Basin-of-attraction analysis, ranked by basin size.
result = m.basins()
print(result.n_fixed_points, result.cycle_count, result.total)
result.to_records()                # flat list of dicts (pandas-friendly)

# Score a candidate scenario and check self-consistency.
m.impact_balance({"WTRD": "FT", "WSEC": "Rlx", "WECO": "RpdGr"})
```

## Editing a model in place (no re-parse)

Tweak individual cross-impact judgements and re-run immediately:

```python
old = m.set_impact(source=("WTRD", "FT"), target=("WSEC", "Alrt"), value=3)
m.consistent_scenarios()           # recomputed on the edited matrix

# Sweep one impact over several values (one parse, N runs); restores after.
sweep_impact(m, ("WTRD", "FT"), ("WSEC", "Alrt"), [-3, 0, 3, 6])

# Branch a model to try edits against a baseline.
variant = m.copy()
variant.set_impact(("WTRD", "FT"), ("WSEC", "Alrt"), 9)
```

Descriptors and variants may be given by name or by 0-based index. Editing
changes the numeric expert judgements (matrix cells); structural changes
(adding/removing descriptors or variants) resize the matrix and require a
reload.

## Batch over many models

```python
records = run_models("path/to/scw_dir", analysis="basins")   # dir, glob, or list
# one record per file: path, n_descriptors, n_scenarios, n_consistent, ...

from crossimpactbalances import to_dataframe
df = to_dataframe(records)          # requires the [pandas] extra
```

## Backends

The same `Model` API runs on either of two engine backends, selected by
`Model.load(..., backend=...)` or the `CIB_BACKEND` environment variable:

| Backend | Ships source? | How it runs |
|---------|---------------|-------------|
| `juliacall` | **Yes** — develops the Julia source package | In-process via juliacall; best for developing in this repo. |
| `native` | **No** | A compiled `libcib` shared library (machine code) driven via `ctypes`. |
| `auto` (default) | — | Uses the compiled library if one is found, else juliacall. |

The two backends each embed a Julia runtime, so only one can be active per
process — pin it rather than mixing.

## Shipping without source (native backend)

To distribute the engine as a black box — no `.jl` files on the user's
machine — build the shared library and bundle it into the wheel.

1. **Build the library** on each platform you ship (Linux/macOS/Windows ×
   x86_64/arm64), from the repo root:

   ```bash
   julia --project=build build/build_library.jl
   # → build/cib-lib/lib/libcib.{so,dylib,dll} (+ headers, runtime)
   ```

2. **Stage it into the package** so the wheel includes it:

   ```bash
   cp -r build/cib-lib python/crossimpactbalances/_libcib
   ```

   `_native.py` looks here first (then `build/cib-lib/`, then the
   `CIB_NATIVE_LIB` env var). `pyproject.toml` already bundles
   `crossimpactbalances/_libcib/**` into the wheel.

3. **Build a platform wheel** (`pip wheel python/` or `python -m build`). The
   resulting wheel carries the compiled engine and needs no Julia and no
   source. Because the runtime is platform-specific, publish one wheel per
   platform (a manylinux/macos/windows tag each).

At runtime the native backend needs no juliacall and no network. Set
`CIB_NATIVE_THREADS=auto` before first use for multi-threaded exhaustive/basin
analysis.

## Packaging note (juliacall backend)

`juliapkg.json` pins the source engine with a `dev` path to the in-repo Julia
source — correct for developing in this repository. To distribute the
juliacall backend on its own, replace the `dev`/`path` entry with a `url`/`rev`
entry pointing at the repository (note this puts the Julia source on the user's
machine — use the native backend if that is not acceptable).

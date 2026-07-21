# Building the Windows installer (.msi)

This directory packages the [desktop app](../app) into a self-contained
Windows `.msi` that installs and runs **without Julia or the source code** on
the target machine.

The build has two stages:

1. **Compile** — [`create_app.jl`](create_app.jl) uses
   [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) to
   turn `app/` into a standalone, relocatable executable bundle
   (`build/CrossImpactBalances-app/` with `bin\`, `lib\`, `share\`, and the
   [`launcher.cmd`](launcher.cmd) that starts the app on all CPU cores).
2. **Package** — [`installer.wxs`](installer.wxs) uses
   [WiX v5](https://wixtoolset.org/) to wrap that bundle into an MSI with a
   Start-menu shortcut and clean upgrade/uninstall behaviour.

> A Julia application cannot be cross-compiled: the compile stage must run on
> **Windows** to produce a Windows executable. The easiest path is the
> automated GitHub Actions workflow below — no local Windows setup required.

## Automated build (recommended)

The workflow [`.github/workflows/build-installer.yml`](../.github/workflows/build-installer.yml)
runs both stages on a `windows-latest` runner and uploads the `.msi`.

- **On demand:** GitHub → *Actions* → *Build Windows installer* → *Run
  workflow*, enter a version (e.g. `0.2.0`). Download the `.msi` from the run's
  *Artifacts*.
- **On a release:** push a tag like `v0.2.0`. The workflow builds the `.msi`
  and attaches it to that GitHub Release, ready to hand out.

## Manual build (on a Windows machine)

Prerequisites: [Julia 1.10+](https://julialang.org/downloads/), a C compiler on
`PATH` (`gcc` — MinGW; PackageCompiler also honours the `JULIA_CC` env var),
and the [.NET SDK](https://dotnet.microsoft.com/download) (for WiX).

From the repository root, in a shell where `%CD%` is that root:

```bat
:: 1. one-time environment setup (Pkg.develop resolves the engine by path)
julia --project=app   -e "using Pkg; Pkg.develop(path=\".\"); Pkg.instantiate()"
julia --project=build -e "using Pkg; Pkg.instantiate()"

:: 2. compile the standalone app bundle (several minutes)
julia --project=build build\create_app.jl "%CD%\build\CrossImpactBalances-app"

:: 3. build the MSI
dotnet tool install --global wix --version 5.*
wix build build\installer.wxs -arch x64 ^
    -d AppDir="%CD%\build\CrossImpactBalances-app" ^
    -d Version=0.2.0 ^
    -o CrossImpactBalances-0.2.0-x64.msi
```

The result, `CrossImpactBalances-0.2.0-x64.msi`, installs to
`Program Files\Cross-Impact Balances` and adds a **Cross-Impact Balances**
Start-menu entry that launches the app in the browser.

## Notes

- **Size.** The bundle embeds the Julia runtime, so the installer is on the
  order of a few hundred MB — expected for a compiled Julia app.
- **Threads.** The Start-menu shortcut runs `launcher.cmd`, which sets
  `JULIA_NUM_THREADS=auto` so the exhaustive search and basin analysis use
  every core. Running `bin\CrossImpactBalances.exe` directly works too, but
  single-threaded unless that variable is set.
- **UpgradeCode.** The GUID in `installer.wxs` is fixed so future versions
  upgrade in place. Keep it stable across releases.
- **Version.** Use the same `Version` (X.Y.Z) as the release tag; MSI requires
  numeric `major.minor.patch`.

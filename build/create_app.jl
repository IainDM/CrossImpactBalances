# Compile the CIBApp desktop application into a standalone, relocatable
# executable bundle with PackageCompiler.
#
#   julia --project=build build/create_app.jl [OUTDIR]
#
# The result (default: build/CrossImpactBalances-app) contains bin/, lib/ and
# share/ and runs with no Julia installation present. On Windows the launched
# executable is bin/CrossImpactBalances.exe; the MSI wraps this bundle.
#
# This step is OS-native: run it on Windows to produce a Windows bundle. The
# GitHub Actions workflow (.github/workflows/build-installer.yml) does exactly
# that on a windows-latest runner.

using PackageCompiler

const ROOT    = normpath(joinpath(@__DIR__, ".."))
const APP_SRC = joinpath(ROOT, "app")
const OUTDIR  = length(ARGS) >= 1 ? ARGS[1] :
                joinpath(ROOT, "build", "CrossImpactBalances-app")

@info "Building standalone app" app = APP_SRC dest = OUTDIR

create_app(
    APP_SRC,
    OUTDIR;
    executables = ["CrossImpactBalances" => "julia_main"],
    precompile_execution_file = joinpath(@__DIR__, "precompile_app.jl"),
    incremental = false,
    filter_stdlibs = false,
    force = true,
    include_lazy_artifacts = false,
)

# Ship the thread-enabling launcher next to the bundle so the Start-menu
# shortcut can start the app with all CPU cores (JULIA_NUM_THREADS=auto).
let src = joinpath(@__DIR__, "launcher.cmd"),
    dst = joinpath(OUTDIR, "Cross-Impact Balances.cmd")
    isfile(src) && cp(src, dst; force = true)
end

@info "Done" bundle = OUTDIR

# Compile the CIBCApi C surface into a standalone, relocatable **shared
# library** with PackageCompiler. The result contains machine code only — no
# Julia source — so the engine can be shipped without exposing its source.
#
#   julia --project=build build/build_library.jl [OUTDIR]
#
# Output (default: build/cib-lib):
#   lib/libcib.{so,dylib,dll}   the compiled engine (+ bundled Julia runtime)
#   include/cib.h, julia_init.h C headers (init_julia / shutdown_julia + cib_*)
#   share/, lib/julia/          runtime support files
#
# The bundle is self-contained and runs with no Julia installation present.
# This step is OS/arch-native: run it on each platform you want to ship
# (Linux/macOS/Windows × x86_64/arm64) to produce that platform's library.

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CAPI = joinpath(ROOT, "capi")
const OUT  = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "build", "cib-lib")

# 1. Make the CIBCApi project resolve CrossImpactBalances from this repo.
@info "Resolving CIBCApi dependencies" capi = CAPI
Pkg.activate(CAPI)
Pkg.develop(path = ROOT)
Pkg.instantiate()

# 2. Switch to the build project (carries PackageCompiler) and compile.
Pkg.activate(joinpath(ROOT, "build"))
Pkg.instantiate()
using PackageCompiler

@info "Building shared library" dest = OUT
create_library(
    CAPI,
    OUT;
    lib_name = "cib",
    precompile_execution_file = joinpath(CAPI, "precompile.jl"),
    incremental = false,
    filter_stdlibs = false,
    force = true,
    include_lazy_artifacts = false,
)

@info "Done" bundle = OUT

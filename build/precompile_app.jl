# Warm-up workload run during compilation so PackageCompiler bakes the hot
# code paths (parsing, exhaustive search, basin analysis, and the JSON/CSV
# response builders) into the app image. It must run without error.

using CrossImpactBalances
using CIBApp

const SAMPLES = normpath(joinpath(@__DIR__, "..", "test", "sample_files"))

for f in ("CIB_global.scw", "bench_typical.scw")
    path = joinpath(SAMPLES, f)
    isfile(path) || continue
    cib = load_scw(path; kernel = Vector{Vector{Int}}())
    find_consistent(cib; exhaustive = true)          # sweep (auto for small spaces)
    # Force the branch-and-bound strategy too: spaces ≥ 10^5 scenarios take it
    # via :auto, and without this line it would be JIT-compiled on the user's
    # first click instead of baked into the image.
    find_consistent(cib; exhaustive = true, algorithm = :bnb)
    find_basins(cib)                                  # two-phase basins
end

# Exercise the app's own request-handling helpers.
let path = joinpath(SAMPLES, "CIB_global.scw")
    if isfile(path)
        cib = load_scw(path; kernel = Vector{Vector{Int}}())
        state = CIBApp.AppState()
        CIBApp.to_json(CIBApp.analyze_consistent(cib))
        CIBApp.to_json(CIBApp.analyze_basins(cib, state))
        CIBApp.load_from_text(read(path, String))
    end
end

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
    find_consistent(cib)                              # whatever :auto picks — the user's path
    # Force BOTH search strategies as well, so neither is left to be JIT-compiled
    # on the user's first click. :auto only picks one of them for a given model
    # size, so relying on the call above would bake in only that one.
    find_consistent(cib; algorithm = :sweep)
    find_consistent(cib; algorithm = :bnb)
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

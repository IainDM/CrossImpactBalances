# Execution trace for PackageCompiler: exercise the engine (the expensive
# code to compile) and the C wrappers so their compiled instances are baked
# into the shared library. Runs in the CIBCApi project context.

using CrossImpactBalances
using CIBCApi

const SAMPLE = joinpath(@__DIR__, "..", "test", "sample_files", "CIB_global.scw")

# --- Engine hot paths ---
cib = load_scw(SAMPLE; kernel = Vector{Vector{Int}}())
find_consistent(cib)                      # whatever :auto picks
find_consistent(cib; algorithm = :sweep)  # bake both search strategies, so
find_consistent(cib; algorithm = :bnb)    # neither is JIT-compiled on first use
find_basins(cib)
u = inv_signature(cib, 0)
impact_balance(cib, u)
succession_step(GlobalSuccession(), cib, u)
signature(cib, u)
set_impact!(cib, 0, 0, 1, 0, 1)
get_impact(cib, 0, 0, 1, 0)

# --- C wrappers (bake the JSON marshalling method instances) ---
GC.@preserve SAMPLE begin
    reply = CIBCApi.cib_load(Cstring(pointer(SAMPLE)))
    CIBCApi.cib_free_string(reply)
end

opts = "{}"
GC.@preserve opts begin
    CIBCApi.cib_free_string(CIBCApi.cib_consistent(Cint(1), Cstring(pointer(opts))))
    CIBCApi.cib_free_string(CIBCApi.cib_basins(Cint(1), Cstring(pointer(opts))))
end
scen = "[0,0,0]"
GC.@preserve scen begin
    CIBCApi.cib_free_string(CIBCApi.cib_impact_balance(Cint(1), Cstring(pointer(scen))))
    CIBCApi.cib_free_string(CIBCApi.cib_signature(Cint(1), Cstring(pointer(scen))))
end
CIBCApi.cib_free(Cint(1))

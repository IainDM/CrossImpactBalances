"""
    CIB

Cross-impact balance analysis object.

Fields:
- `descriptors`: ordered list of descriptor names
- `variants`: dict mapping descriptor name to list of variant names
- `numberOfVariants`: number of variants per descriptor (can be different for each descriptor)
- `cim`: the cross-impact matrix (n × n), n = sum of all variants. Row `r` holds the impacts contributed by variant `r` onto every other variant.
- `cim_t`: the transpose of `cim` — `cim_t[j, r] == cim[r, j]`. Stored alongside `cim` so that "sum a fixed-row vector across many `j`" hot loops (impact_balance, find_basins) can iterate contiguously through column-major storage, which unlocks SIMD. Construction cost is one `permutedims` at load time.
- `numberOfDimensions`: total number of variants (size of CIM)
- `numberOfDescriptors`: number of descriptors
- `consistentScenarios`: list of consistent scenarios (each a Vector{Int}, 0-based variant indices)
- `desc_offsets`: cumulative 0-based variant offsets per descriptor, shows where in the CIM each descriptor's variants start

A `struct` in Julia is immutable: the field *bindings* can never be reassigned after construction (like a C# readonly record or a frozen dataclass). The arrays the fields point at can still have their *contents*
changed — which is how `consistentScenarios` gets filled in after the object is built.
"""
struct CIB
    descriptors::Vector{String}
    variants::Dict{String, Vector{String}}
    numberOfVariants::Vector{Int}
    cim::Matrix{Int}
    cim_t::Matrix{Int}
    numberOfDimensions::Int
    numberOfDescriptors::Int
    consistentScenarios::Vector{Vector{Int}}
    desc_offsets::Vector{Int}
end

"""
    _score_type(cib) -> Type{<:Signed}

Narrowest safe integer type for impact-balance accumulation. `Int16` when `4 * (ndesc + 1) * maximum(abs, cim)` fits (margin covers the incremental row-delta updates; Int16 arithmetic wraps silently, so this guard is the only protection), otherwise `Int`. Realistic CIB matrices (entries ±3) are far inside the Int16 bound. Narrower integers matter because SIMD processes twice as many Int16 values per instruction as Int32.
    
Basically, Int16 is faster so use it if possible, and only for extremely large CIB matrices would the value of a scenario possibly be above the Int16 limit of 32,767)
"""
function _score_type(cib::CIB)
    # `cond ? a : b` is the ternary operator, as in C#.
    maxAbsoluteImpact = cib.numberOfDimensions == 0 ? 0 : Int(maximum(abs, cib.cim))
    return 4 * (cib.numberOfDescriptors + 1) * maxAbsoluteImpact <= Int(typemax(Int16)) ? Int16 : Int
end
#TODO: should really be a property of the CIB not computed multiple times. also, it's in the wrong place, this section is for sweeps

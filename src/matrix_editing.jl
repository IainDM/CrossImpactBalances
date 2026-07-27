# Everything else in this file addresses the matrix by flat row/column number,
# because that is what the hot loops need. These helpers are the human-facing
# way in: they let a caller say "the impact of Trade=Free on Growth=High"
# instead of "cim[4, 9]", and they are what the Python and C wrappers call.
#
# Nothing inside the engine uses them — they exist purely so that a model can
# be inspected and edited from outside without re-parsing the .scw file, which
# is what makes sensitivity sweeps (change one judgement, re-run) cheap.

"""
    _desc_index(cib, descriptor) -> Int

Resolve a descriptor given by name (`AbstractString`) or 0-based index (`Integer`) to its 1-based position in `cib.descriptors`. Throws on an unknown name or an out-of-range index.
"""
# Two methods, one name: Julia picks between them on the TYPE of the second
# argument (multiple dispatch), so the caller can pass either "Trade" or 0.
function _desc_index(cib::CIB, descriptor::AbstractString)
    position = findfirst(==(String(descriptor)), cib.descriptors)
    isnothing(position) && throw(ArgumentError(
        "Unknown descriptor: \"$descriptor\". Available: $(join(cib.descriptors, ", "))"))
    return position
end

function _desc_index(cib::CIB, descriptor::Integer)
    (0 <= descriptor < cib.numberOfDescriptors) || throw(ArgumentError(
        "Descriptor index $descriptor out of range 0:$(cib.numberOfDescriptors - 1)"))
    return Int(descriptor) + 1     # callers count descriptors from 0, arrays from 1
end

"""
    _table_index(cib, descriptor, variant) -> Int

Resolve a (descriptor, variant) pair to the 1-based flat row/column index into `cib.cim`. Both may be given by name (`AbstractString`) or by 0-based index (`Integer`).

This is the same offset arithmetic the scoring loops do inline — descriptor `i`'s variants start at `desc_offsets[i]`, so variant `v` of that descriptor lives at `desc_offsets[i] + v + 1`.
"""
function _table_index(cib::CIB, descriptor, variant::AbstractString)
    descriptorPosition = _desc_index(cib, descriptor)
    descriptorName = cib.descriptors[descriptorPosition]
    variantNames = cib.variants[descriptorName]
    variantPosition = findfirst(==(String(variant)), variantNames)
    isnothing(variantPosition) && throw(ArgumentError(
        "Unknown variant \"$variant\" for descriptor \"$descriptorName\". " *
        "Available: $(join(variantNames, ", "))"))
    # findfirst already returns a 1-based position, so no +1 here.
    return cib.desc_offsets[descriptorPosition] + variantPosition
end

function _table_index(cib::CIB, descriptor, variant::Integer)
    descriptorPosition = _desc_index(cib, descriptor)
    variantCount = cib.numberOfVariants[descriptorPosition]
    (0 <= variant < variantCount) || throw(ArgumentError(
        "Variant index $variant out of range 0:$(variantCount - 1) " *
        "for descriptor \"$(cib.descriptors[descriptorPosition])\""))
    return cib.desc_offsets[descriptorPosition] + Int(variant) + 1
end

"""
    set_impact!(cib, src_desc, src_var, tgt_desc, tgt_var, value) -> Int

Set the cross-impact contributed by the source variant (`src_desc`=`src_var`) onto the target variant (`tgt_desc`=`tgt_var`), i.e. `cim[source, target]`, and return the previous value.

Descriptors and variants may be given by name (`AbstractString`) or by 0-based index (`Integer`). This updates **both** the cross-impact matrix `cim` and its stored transpose `cim_t`, so that every analysis routine — which may read either — sees a consistent matrix.

Because the matrices are edited in place, the already-loaded model changes without any re-parse of the `.scw` file; the next [`find_consistent`](@ref) or [`find_basins`](@ref) call reflects the new value. That is what makes it cheap to sweep one expert judgement across a range and watch the consistent scenarios move.

Note the `!` in the name: by Julia convention it warns that the argument is modified rather than copied.
"""
function set_impact!(cib::CIB, src_desc, src_var, tgt_desc, tgt_var, value::Integer)
    sourceIndex = _table_index(cib, src_desc, src_var)
    targetIndex = _table_index(cib, tgt_desc, tgt_var)
    previousValue = cib.cim[sourceIndex, targetIndex]
    newValue = Int(value)
    cib.cim[sourceIndex, targetIndex]   = newValue
    cib.cim_t[targetIndex, sourceIndex] = newValue   # keep the transpose in sync
    return previousValue
end

"""
    get_impact(cib, src_desc, src_var, tgt_desc, tgt_var) -> Int

Return the current cross-impact `cim[source, target]` contributed by the source variant (`src_desc`=`src_var`) onto the target variant (`tgt_desc`=`tgt_var`). Descriptors and variants may be given by name or 0-based index. The read-only partner of [`set_impact!`](@ref).
"""
function get_impact(cib::CIB, src_desc, src_var, tgt_desc, tgt_var)
    sourceIndex = _table_index(cib, src_desc, src_var)
    targetIndex = _table_index(cib, tgt_desc, tgt_var)
    return cib.cim[sourceIndex, targetIndex]
end

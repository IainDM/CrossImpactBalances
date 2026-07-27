"""
    impact_balance(cib, u) -> Vector{Int}

Compute the impact balance vector for scenario `u`: one score per variant across all descriptors (length `numberOfDimensions`). The score of a variant is the total impact the scenario's chosen variants exert on it — the higher the score, the better that variant fits the scenario.
"""
function impact_balance(cib::CIB, scenario::Vector{Int})
    numberOfDimensions = cib.numberOfDimensions
    impactBalance = zeros(Int, numberOfDimensions)

    # Work with the transposed matrix. this is a trick to allow the CPU to use SIMD instructions.
    cimTranspose = cib.cim_t
    descriptorOffsets = cib.desc_offsets

    # @inbounds tells Julia to skip array bounds checking inside the block — the equivalent of unchecked array access in C#. Safe here because all the indices are derived from the CIB's own dimensions.
    @inbounds for (descriptorIndex, chosenVariant) in enumerate(scenario)
        # The matrix row belonging to this descriptor's chosen variant (+1 converts the 0-based variant number to a 1-based array index).
        sourceRow = descriptorOffsets[descriptorIndex] + chosenVariant + 1
        # @simd allows the compiler to vectorise this accumulation loop (process several array elements per CPU instruction).
        @simd for targetVariant in 1:numberOfDimensions
            impactBalance[targetVariant] += cimTranspose[targetVariant, sourceRow]
        end
    end
    return impactBalance
end

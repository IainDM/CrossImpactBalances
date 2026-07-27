"""
    signature(cib, u) -> Int

Compute a unique, sequential, integer signature for scenario `u` (0-based variant indices).

A scenario is a choice of one variant per descriptor. Treating those choices as the digits of a mixed-radix number (each descriptor's "digit" can count up to its own variant count) maps every scenario to a distinct integer in `0:max_signature(cib)` — exactly like reading [7, 2, 4] as a decimal number, except each position can have a different base. This lets the search code step through all scenarios with a single counter and avoids expensive modulus operations
"""
function signature(cib::CIB, scenario::Vector{Int})
    signatureValue = 0
    placeValue = 1     # the "place value" of the current digit (1, then ×base each step)
    # zip pairs each descriptor's chosen variant with its variant count, like Python's zip.
    for (chosenVariant, variantCount) in zip(scenario, cib.numberOfVariants)
        signatureValue += placeValue * chosenVariant
        placeValue *= variantCount
    end
    return signatureValue
end

"""
    inv_signature(cib, s) -> Vector{Int}

Convert a signature back to a scenario (0-based variant indices), the reverse of the function above
"""
function inv_signature(cib::CIB, signatureValue::Int)
    scenario = Int[]
    for variantCount in cib.numberOfVariants
        # `%` is remainder, `÷` is integer division (Python's //, C#'s /).
        push!(scenario, signatureValue % variantCount)
        signatureValue = signatureValue ÷ variantCount
    end
    return scenario
end

"""
    max_signature(cib) -> Int

The maximum signature value (= total scenarios - 1).
"""
function max_signature(cib::CIB)
    # The largest signature belongs to the scenario picking the last variant
    # of every descriptor. `[... for ... in ...]` is a comprehension, as in
    # Python.
    return signature(cib, [variantCount - 1 for variantCount in cib.numberOfVariants])
end

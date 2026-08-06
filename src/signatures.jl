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
    scenario_count(cib) -> Int128

The total number of scenarios in the model — the product of every descriptor's variant count.

Returned as an `Int128` deliberately. The product of even modest variant counts outruns `Int64` at around forty ternary descriptors, and worse, `Int64` arithmetic *wraps silently* when it overflows — a model would quietly report a garbage size. (It also outruns `Float64` exactness far earlier, at 2⁵³ ≈ 9×10¹⁵: real matrices exist whose scenario count displayed through a double comes out with its final digits rounded.) `Int128` is exact out to ~1.7×10³⁸, which no conceivable CIB matrix approaches, and this function guards even that.

This is the number to display; `max_signature` (= this minus one, as an `Int`) is the number to *iterate to*, and it throws for models whose signatures no longer fit in `Int64` — see its docstring for what still works in that regime.
"""
function scenario_count(cib::CIB)
    totalScenarios = Int128(1)
    for variantCount in cib.numberOfVariants
        # Check BEFORE multiplying: would this multiply pass typemax(Int128)?
        # (`÷` is integer division. Guarding this way needs no wider type.)
        if variantCount > 1 && totalScenarios > typemax(Int128) ÷ variantCount
            throw(ArgumentError(
                "scenario_count: the scenario space exceeds typemax(Int128) ≈ 1.7e38 " *
                "— this model cannot be analysed (or meaningfully counted)"))
        end
        totalScenarios *= variantCount
    end
    return totalScenarios
end

"""
    _signature128(cib, u) -> Int128

The signature of scenario `u` as an `Int128` — the same mixed-radix number [`signature`](@ref) computes, in an integer type wide enough for any model [`scenario_count`](@ref) accepts. The sampling analysis ([`estimate_basins`](@ref)) keys its tallies on this, so it keeps working on models too large for `Int64` signatures, where `signature` itself would silently wrap.
"""
function _signature128(cib::CIB, scenario::Vector{Int})
    signatureValue = Int128(0)
    placeValue = Int128(1)
    for (chosenVariant, variantCount) in zip(scenario, cib.numberOfVariants)
        signatureValue += placeValue * chosenVariant
        placeValue *= variantCount
    end
    return signatureValue
end

"""
    max_signature(cib) -> Int

The maximum signature value (= total scenarios - 1).

Throws an `ArgumentError` for models with more than `typemax(Int64)` ≈ 9.2×10¹⁸ scenarios (for scale: ~130 binary or ~40 ternary descriptors), because signatures beyond that cannot be represented and every signature-driven analysis — [`find_consistent`](@ref), [`find_basins`](@ref) — would silently wrap. Before this guard, `max_signature` itself wrapped silently. [`scenario_count`](@ref) still counts such models exactly, and [`estimate_basins`](@ref) can still estimate their basin shares (it never touches `Int64` signatures); [`influence_structure`](@ref) and [`fix_descriptor`](@ref) can cut them down to exactly analysable pieces.
"""
function max_signature(cib::CIB)
    totalScenarios = scenario_count(cib)
    if totalScenarios - 1 > Int128(typemax(Int))
        throw(ArgumentError(
            "max_signature: this model has $(totalScenarios) scenarios, more than " *
            "typemax(Int64) = $(typemax(Int)) — signatures cannot be represented, so the " *
            "exact signature-based analyses are unavailable. scenario_count(cib) still " *
            "counts the space exactly; estimate_basins(cib) still estimates basin shares; " *
            "influence_structure(cib) / fix_descriptor(cib, d, v) can reduce the model to " *
            "exactly analysable pieces."))
    end
    # Identical value to the old comprehension-and-signature computation
    # (the largest signature belongs to the scenario picking the last variant
    # of every descriptor), minus the silent wraparound past typemax(Int64).
    return Int(totalScenarios - 1)
end

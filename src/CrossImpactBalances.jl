"""
    CrossImpactBalances

Cross-Impact Balance (CIB) scenario analysis. A Julia implementation of the methodology in Weimer-Jehle (2006), ported from the Python [sei-international/cibsa](https://github.com/sei-international/cibsa) library with a two-phase basin-of-attraction analysis, a threaded incremental exhaustive sweep, and an exact branch-and-bound search that prunes provably-inconsistent regions of the scenario space. The succession dynamics can be changed through [`SuccessionRule`](@ref). There is a default called [`GlobalSuccession`](@ref) which follows the standard ScenarioWizard rule.

# Reference
Weimer-Jehle, W. (2006). Cross-impact balances: A system-theoretical approach to cross-impact analysis. *Technological Forecasting and Social Change*, 73(4), 334–361.
"""
module CrossImpactBalances

# `export` lists the names visible after `using CrossImpactBalances` — like
# Python's `__all__` or C#'s `public`. Everything else in the module is still
# reachable as `CrossImpactBalances.name`, it just isn't brought into scope
# automatically. Names starting with an underscore are internal helpers.
export CIB, load_scw, load_solutions,
       impact_balance,
       set_impact!, get_impact,
       SuccessionRule, GlobalSuccession, SequentialSuccession,
       succession_step,
       find_consistent, find_basins,
       signature, inv_signature, max_signature

# ── Julia notes for readers coming from Python or .NET ──────────────────────
#
# Features that recur throughout this file:
#
# * 1-based indexing: Julia arrays start at index 1 (like MATLAB/Fortran),
#   and a range `1:n` includes BOTH endpoints. Scenarios themselves store
#   0-based variant numbers (matching the Python original and the signature
#   arithmetic), which is why many array accesses add `+ 1`.
#
# * Multiple dispatch: several functions below share one name but differ in
#   their argument TYPES (e.g. `succession_step(::GlobalSuccession, ...)` vs
#   `succession_step(::SequentialSuccession, ...)`). Julia picks the method
#   from the runtime types of all arguments — like C# method overloading,
#   but resolved dynamically, and it is the standard way to write what
#   Python/C# would express with an interface or virtual method.
#
# * A trailing `!` in a function name (e.g. `push!`, `_successor_table!`) is
#   a naming convention meaning "this function mutates one of its
#   arguments". It has no effect on behaviour — it is a warning label.
#
# * Words starting with `@` are macros — code transformations applied before
#   compilation. The ones used here are performance/threading annotations:
#   `@inbounds` (skip array bounds checks), `@simd` (allow vectorisation),
#   `Threads.@spawn` (run on a worker thread, like Task.Run), `@sync` (wait
#   for all tasks spawned inside the block, like Task.WaitAll), and `@view`
#   (slice an array without copying it).
#
# * `:auto`, `:sweep`, `:bnb` are Symbols — lightweight interned names, used
#   here the way C# would use an enum or Python a short string constant.
#
# * `Union{Nothing, Int}` is a nullable type (C#'s `int?`); `nothing` plays
#   the role of null/None, tested with `isnothing(x)` or `x === nothing`.
#
# * `condition && action` and `condition || action` are used as one-line if
#   statements: `&&` runs the action only when the condition is true,
#   `||` only when it is false (e.g. `isempty(s) && continue`).
#
# * `"""docstrings"""` directly above a definition are attached to it as
#   documentation (visible via `?name` in the REPL); `[`Name`](@ref)` inside
#   them is a cross-reference link for the documentation generator.
#
# * `#region` / `#endregion` comments are only editor folding markers
#   (recognised by VS Code); Julia ignores them.

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



#region "files parsers"
# ─── .scw file parser ───────────────────────────────────────────────────────

"""
    load_scw(scw_file; sl_file=nothing, kernel=nothing, algorithm=:auto) -> CIB

Parse a ScenarioWizard .scw file and optionally a .sl solutions file.
Returns a fully populated CIB object.

Unless a `kernel` or an `sl_file` is supplied, the kernel is computed by [`find_consistent`](@ref), which always searches the full scenario space and uses every available thread (start Julia with `julia -t auto`). `algorithm` is forwarded to select the search strategy.
"""
function load_scw(scw_file::String; sl_file::Union{String,Nothing}=nothing,
                  kernel::Union{Vector{Vector{Int}},Nothing}=nothing,
                  algorithm::Symbol=:auto)
    descriptors = String[]
    variants = Dict{String, Vector{String}}()
    variantCounts = Int[]        # how many variants each descriptor has

    # The .scw format is line-oriented: a header section listing descriptors
    # ('&' lines) and their variants ('-' lines), then six '#'-separated
    # sections, the sixth of which is the cross-impact matrix as CSV rows.
    # We read it with a small state machine: parserState 0 is the header,
    # each '#' line advances the state, and state 5 collects matrix rows.
    parserState = 0
    totalVariants = 0            # running count of all variants seen
    descriptorsSeen = -1         # -1 until the first descriptor line arrives
    variantsInCurrent = 0        # variants seen for the descriptor being read
    currentDescriptor = ""

    matrixRows = Vector{Vector{Int}}()   # raw CIM rows, validated below

    for line in eachline(scw_file)
        stripped = lstrip(line)
        isempty(stripped) && continue    # skip blank lines

        if parserState == 0
            if stripped[1] == '&'
                # A new descriptor. Before starting it, record how many
                # variants the previous descriptor had (skipped for the very
                # first one, when descriptorsSeen is still -1).
                descriptorName = strip(stripped[2:end])
                push!(descriptors, descriptorName)
                variants[descriptorName] = String[]
                if descriptorsSeen > -1
                    push!(variantCounts, variantsInCurrent)
                end
                variantsInCurrent = 0
                descriptorsSeen += 1
                currentDescriptor = descriptorName
            elseif stripped[1] == '-'
                # A variant belonging to the descriptor currently being read.
                push!(variants[currentDescriptor], strip(stripped[2:end]))
                totalVariants += 1
                variantsInCurrent += 1
            elseif stripped[1] == '#'
                # End of the header: close out the last descriptor's count
                # and move to section 1.
                parserState = 1
                push!(variantCounts, variantsInCurrent)
            end
        elseif parserState < 5 && stripped[1] == '#'
            # Sections 1-4 hold ScenarioWizard metadata we don't need; just
            # count the '#' separators until the matrix section arrives.
            parserState += 1
        elseif parserState == 5
            if stripped[1] == '#'
                parserState += 1     # '#' after the matrix: we're done reading
            else
                # A matrix row: comma-separated integers. `parse.(Int, ...)`
                # uses broadcasting — the dot applies `parse` to every element
                # of the split list at once (like Python's map / LINQ Select).
                rowValues = parse.(Int, split(stripped, ','))
                push!(matrixRows, rowValues)
            end
        end
    end

    # `$(...)` inside a string is interpolation, like Python f-strings.
    totalVariants == 0 && error("load_scw: no variants found in $(scw_file) — file is empty or malformed")
    length(matrixRows) == totalVariants ||
        error("load_scw: cross-impact matrix has $(length(matrixRows)) " *
              "rows but $totalVariants variants in $(scw_file)")

    # Copy the validated rows into a proper square matrix.
    cim = zeros(Int, totalVariants, totalVariants)
    for (rowIndex, rowValues) in enumerate(matrixRows)
        length(rowValues) == totalVariants ||
            error("load_scw: row $rowIndex of CIM has $(length(rowValues)) entries but $totalVariants expected")
        for (columnIndex, value) in enumerate(rowValues)
            cim[rowIndex, columnIndex] = value
        end
    end

    numberOfDescriptors = descriptorsSeen + 1

    # desc_offsets[i] is where descriptor i's block of variants starts in the
    # matrix (0-based). E.g. with variant counts [3, 2, 4] the offsets are
    # [0, 3, 5]: descriptor 2's variants occupy matrix rows 4 and 5.
    desc_offsets = Vector{Int}(undef, numberOfDescriptors)  # undef = allocate without initialising
    runningOffset = 0
    for descriptorIndex in 1:numberOfDescriptors
        desc_offsets[descriptorIndex] = runningOffset
        runningOffset += variantCounts[descriptorIndex]
    end

    # Precompute the transpose so impact_balance / find_basins can do
    # contiguous SIMD column reads in cim_t (= the row vectors of cim).
    # Julia stores matrices column-major (like Fortran, unlike C#/NumPy's
    # default row-major), so summing down a column is the fast direction.
    cim_t = permutedims(cim)

    # Build the CIB with an empty kernel first, then fill the kernel in
    # place (the struct is immutable, but the vector's contents are not).
    cib = CIB(descriptors, variants, variantCounts, cim, cim_t,
              totalVariants, numberOfDescriptors,
              Vector{Vector{Int}}(), desc_offsets)

    if !isnothing(kernel)
        append!(cib.consistentScenarios, kernel)
    elseif !isnothing(sl_file)
        append!(cib.consistentScenarios, load_solutions(cib, sl_file))
    else
        append!(cib.consistentScenarios, find_consistent(cib; algorithm=algorithm))
    end

    return cib
end

# ─── .sl file parser ────────────────────────────────────────────────────────

"""
    load_solutions(cib::CIB, sl_file::String) -> Vector{Vector{Int}}

Parse a ScenarioWizard .sl solutions file. Returns 0-based variant indices.
"""
function load_solutions(cib::CIB, sl_file::String)
    solutions = Vector{Vector{Int}}()
    for line in eachline(sl_file)
        stripped = lstrip(line)
        isempty(stripped) && continue
        stripped[1] != '"' && continue   # solution lines start with a quote

        # Extract the quoted index string, e.g. "2 3 2"
        quoted = match(r"^\"([^\"]+)\"", stripped)
        isnothing(quoted) && continue

        indices = parse.(Int, split(strip(quoted.captures[1])))
        # ScenarioWizard numbers variants from 1; internally we use 0-based
        # numbers. `.- 1` is a broadcast: subtract 1 from every element.
        push!(solutions, indices .- 1)
    end
    return solutions
end
#endregion

#region "─── Signature (unique integer ID for a scenario) ──────────────────────────"

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
#endregion

#region "Calculate Impact balance"

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

#endregion

#region "Reading and editing the cross-impact matrix"

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

#endregion

#region ─── Succession

"""
    SuccessionRule

Succession is one of the two core concepts of CIB - "how do we move from one scenario to the next".
The other core concept is consistency - "is this scenario better than others".
Note that these concepts are separate in that the question of whether there is a better scenario is separate from the question of how to get there.

To allow flexibility in how succession is calculated, we define here an abstract supertype called [`SuccessionRule`](@ref). An abstract type cannot be instantiated itself — it only serves as a parent for concrete rule types, the way an interface or abstract base class would in C# or Python.

To add in a new rule, define

    struct MyRule <: SuccessionRule end

(`<:` means "is a subtype of") and a single method

    succession_step(rule::MyRule, cib::CIB, u::Vector{Int}) -> Vector{Int}

and every analysis routine — [`find_consistent`](@ref) and [`find_basins`](@ref) — works with it immediately, because they dispatch on the rule's type. The only assumption is that a SuccessionRule depends only on the current scenario, not the path taken to get there (in technical terms, it's a Markov process).
"""
abstract type SuccessionRule end

"""
The standard succession approach: in one step, every descriptor is moved to the variant with the highest impact score given the *current* scenario. Ties favour the current variant, then the lowest index.
"""
struct GlobalSuccession <: SuccessionRule end

"""
Sequential (successive / Gauss–Seidel) succession: descriptors are updated one at a time in descriptor order, each using the impact balance of the scenario *as already partially updated within the same step*.
"""
struct SequentialSuccession <: SuccessionRule end

"""
    succession_step(cib, u) -> Vector{Int}
    succession_step(rule, cib, u) -> Vector{Int}

One step of succession under `rule`.
"""
function succession_step(::GlobalSuccession, cib::CIB, scenario::Vector{Int})

    # Score every variant against the CURRENT scenario once, then let each descriptor independently pick its best variant from those scores.
    impactBalance = impact_balance(cib, scenario)
    successor = copy(scenario)
    firstColumn = 1     # start of this descriptor's block in impactBalance (1-based). will be increased by the number of variants to point to the start of the next descriptor's block later on

    for descriptorIndex in 1:cib.numberOfDescriptors
        variantCount = cib.numberOfVariants[descriptorIndex]
        lastColumn = firstColumn + variantCount - 1

        # @view slices without copying — like a NumPy view or a C# Span, not the fresh array a plain slice would allocate.
        descriptorScores = @view impactBalance[firstColumn:lastColumn]

        # Seed the running best with the CURRENT variant's score, so a tie  can never displace it.
        bestScore = descriptorScores[scenario[descriptorIndex] + 1]
        for variantIndex in 0:variantCount-1

            # Strict > means ties keep the current variant, and among equal challengers the lowest index wins (it is reached first).
            if descriptorScores[variantIndex + 1] > bestScore
                bestScore = descriptorScores[variantIndex + 1]
                successor[descriptorIndex] = variantIndex
            end
        end
        firstColumn = lastColumn + 1
    end
    return successor
end

"""
as above except for the SequentialSuccession rule
"""
function succession_step(::SequentialSuccession, cib::CIB, scenario::Vector{Int})
    successor = copy(scenario)
    @inbounds for descriptorIndex in 1:cib.numberOfDescriptors
        # Unlike the global rule, the scores are recomputed for EVERY descriptor, from the partially-updated successor — so descriptor 2 already sees whatever descriptor 1 just changed.
        impactBalance = impact_balance(cib, successor)
        offset = cib.desc_offsets[descriptorIndex]
        variantCount = cib.numberOfVariants[descriptorIndex]

        # Current variant's score seeds the best, so ties favour it.
        bestScore = impactBalance[offset + successor[descriptorIndex] + 1]
        for variantIndex in 0:variantCount-1
            if impactBalance[offset + variantIndex + 1] > bestScore
                bestScore = impactBalance[offset + variantIndex + 1]
                successor[descriptorIndex] = variantIndex
            end
        end
    end
    return successor
end

# Convenience method: calling succession_step without naming a rule means the
# standard global rule. (This one-line `f(x) = ...` form is just compact
# syntax for defining a function — the body is the expression after the `=`.)
succession_step(cib::CIB, scenario::Vector{Int}) =
    succession_step(GlobalSuccession(), cib, scenario)

#endregion

#region "─── Find consistent scenarios ──────────────────────────────────────────────"

"""
    fixed_point_margin(rule) -> Union{Nothing,Int}

A standard fixed point (consistent scenario) is one where no descriptor can be improved by choosing a different variant. By declaring a `fixed_point_margin`, a rule changes this to "cannot be improved *by more than the margin*". This reflects the fact that social systems have some inertia or resistance to change, and so may not shift for a relatively small improvement. `fixed_point_margin = 0` gives the standard behaviour.

Declaring a margin is what lets a rule use the fast search strategies in [`find_consistent`](@ref); a rule that returns `nothing` (the default) is searched by the generic one-scenario-at-a-time scan instead.
"""
fixed_point_margin(::SuccessionRule) = nothing
fixed_point_margin(::GlobalSuccession) = 0
# Sequential succession's fixed points coincide with global's (can be proven by induction) so it can also use the margin-0 fast searches.
fixed_point_margin(::SequentialSuccession) = 0

"""
    find_consistent(cib; rule=GlobalSuccession(), algorithm=:auto,
                    bnb_node_budget=nothing) -> Vector{Vector{Int}}

Find every consistent scenario — every fixed point of the succession map.

`rule` selects the succession *dynamics* ([`SuccessionRule`](@ref); default [`GlobalSuccession`](@ref)).

`algorithm` selects the search strategy (only for rules with a [`fixed_point_margin`](@ref)):
- `:sweep` — enumerate every scenario with the incremental odometer sweep.
- `:bnb`   — branch-and-bound: assign descriptors depth-first and prune subtrees that provably contain no fixed point. Exact — returns the identical kernel — and typically visits a small fraction of the space.
- `:auto`  (default) — the sweep for small spaces (< 10^5 scenarios), otherwise branch-and-bound with a node budget of `n ÷ 16`; if pruning is too weak to pay off, the budget trips and the sweep runs instead.

The returned kernel is ordered by ascending signature.
"""
function find_consistent(cib::CIB; rule::SuccessionRule=GlobalSuccession(),
                         algorithm::Symbol=:auto,
                         bnb_node_budget::Union{Nothing,Int}=nothing)

    return _find_kernel(rule, cib; algorithm=algorithm,
                              bnb_node_budget=bnb_node_budget)
end

"""
    _find_kernel(rule, cib; algorithm, bnb_node_budget)

Exhaustive fixed-point search under `rule`. A rule that declares a [`fixed_point_margin`](@ref) gets the fast threaded sweep / branch-and-bound — the margin parameterises the per-descriptor fixed-point test, so both paths serve global succession (`m = 0`), inertial/threshold rules (`m > 0`) and any other rule with the same separable structure. A rule with no margin falls back to a generic ascending-signature scan that tests `succession_step(rule, u) == u`.
This is because if we know that we can use the standard 'no variant improves the scenario by more than a fixed amount' rule to identify fixed points, there are shortcuts we can use. 
The generic path is left for other rules which might define fixed points differently.
"""
function _find_kernel(rule::SuccessionRule, cib::CIB; algorithm::Symbol=:auto,
                            bnb_node_budget::Union{Nothing,Int}=nothing)
    margin = fixed_point_margin(rule)

    # Generic fallback for rules without a declared margin: walk every scenario in ascending signature order and keep the ones that are their own successor.
    # Single-threaded, one succession_step call per scenario. This is a slow but guaranteed correct for any succession rule approach.
    # The built in rules never come here because both declare margin = 0 so can use the fast path    
    if margin === nothing
        algorithm === :auto || throw(ArgumentError(
            "algorithm=$(repr(algorithm)) needs a rule with a fixed_point_margin; " *
            "this rule uses the generic scan (algorithm=:auto)"))
        consistentScenarios = Vector{Vector{Int}}()
        for currentSignature in 0:max_signature(cib)
            scenario = inv_signature(cib, currentSignature)
            #if this scenario is its own successor then add it to the stack of consistent scenarios
            if succession_step(rule, cib, scenario) == scenario
                push!(consistentScenarios, scenario)
            end
        end
        return consistentScenarios
    end

    # at this point we know the rule declared a margin so we can use the fast searches
    # there's a fast search that checks all scenarios and one that uses 'branch and bound' to only check a subset
    # first check that the algorithm is a valid one
    algorithm in (:auto, :bnb, :sweep) ||
        throw(ArgumentError("algorithm must be :auto, :bnb or :sweep, got $(repr(algorithm))"))
    numberOfScenarios = max_signature(cib) + 1

    # where the caller insists on a full sweep, or the number of scenarios is too low to make it worthwhile
    # doing a branch and bound, just run through them all in order
    if algorithm == :sweep || (algorithm == :auto && numberOfScenarios < 100_000)
        # Note the sweep depends only on the margin, not on the rule itself.
        return _find_kernel_checkall_fast(cib; margin=margin)
    end

    # ok, now we're doing branch and bound to reduce the number of scenarios we need to review
    # this works by gradually introducing each descriptor, and then seeing if it can prove that even if it added every other descriptor, nothing could be a fixed point
    suffixMin, suffixMax = _bnb_bounds(cib)

    # `something(a, b)` returns the first argument that isn't `nothing` — like C#'s ?? null-coalescing operator.
    nodeBudget = something(bnb_node_budget,
                           algorithm == :bnb ? typemax(Int) : numberOfScenarios ÷ 16)
    
    # Get the set of fixed points / consistent scenarios.
    kernel, _ = _bnb_fixed_points(cib, suffixMin, suffixMax;
                                  node_budget=nodeBudget, margin=margin)

    # A `nothing` kernel does NOT mean "no consistent scenarios" — a model with
    # none of those returns an empty list. It means branch-and-bound gave up:
    # it burned through its node budget without finishing, because on this
    # particular matrix the pruning was too weak to be worth it (see the note
    # on the budget in the branch-and-bound section below).
    if !isnothing(kernel)
        return kernel
    end

    # Branch-and-bound gave up, so fall back to checking every scenario.
    return _find_kernel_checkall_fast(cib; margin=margin)
end

#endregion

#region "Exhaustive review of scenarios looking for fixed points"

"""
    _find_kernel_checkall_fast(cib) -> Vector{Vector{Int}}

Splits the scenario space into contiguous chunks scheduled as tasks. Within a chunk, the impact-balance vector is maintained *incrementally*: a mixed-radix odometer advances the scenario, and each digit change updates the balance by the difference of two CIM rows (a `@simd` loop over the narrowed transpose from [`_score_type`](@ref)), so no per-scenario score recomputation happens at all.
"""
function _find_kernel_checkall_fast(cib::CIB; margin::Int=0)
    scoreType = _score_type(cib)
    # `Matrix{scoreType}(...)` copies cim_t into a matrix of the narrower
    # element type, so the hot loops move less memory.
    return _sweep_fixed_points_all(cib, Matrix{scoreType}(cib.cim_t); margin=margin)
end

# `where {ScoreInt<:Signed}` declares a type parameter — like a C# generic method `Sweep<T>(...) where T : ...`.
# Julia compiles a specialised version of the function for each concrete element type it is called with (Int16 or Int here).
# Scores can be negative, so the type must be a signed integer.
function _sweep_fixed_points_all(cib::CIB, cimTranspose::Matrix{ScoreInt};
                             margin::Int=0) where {ScoreInt<:Signed}
    numberOfScenarios = max_signature(cib) + 1

    # divide up all the scenarios into 'chunks' and process each one individually
    # we need to be careful dividing up into threads since the number of scenarios might not be an exact multiple of the number of threads
    numberOfChunks = max(1, min(numberOfScenarios, 16 * Threads.nthreads()))
    chunkSize = cld(numberOfScenarios, numberOfChunks)
    numberOfChunks = cld(numberOfScenarios, chunkSize)

    # One private output vector per chunk, so the threads never write to shared state and need no locks.
    chunkResults = [Vector{Vector{Int}}() for _ in 1:numberOfChunks]
  
    # @sync waits for every task spawned inside the block to finish (like Task.WaitAll); Threads.@spawn schedules the call on a worker thread (like Task.Run).
    @sync for chunkIndex in 1:numberOfChunks
        chunkOutput = chunkResults[chunkIndex]
        firstSignature = (chunkIndex - 1) * chunkSize
        lastSignature = min(chunkIndex * chunkSize, numberOfScenarios) - 1
        Threads.@spawn _sweep_chunk_all!(chunkOutput, cimTranspose,
                                     firstSignature, lastSignature,
                                     cib.numberOfVariants, cib.desc_offsets,
                                     cib.numberOfDescriptors, cib.numberOfDimensions,
                                     margin)
    end

    # Chunks are in scenario number order so we can just concatenate them, no need to dedupe or to sort
    kernel = Vector{Vector{Int}}()
    for chunkIndex in 1:numberOfChunks
        append!(kernel, chunkResults[chunkIndex])
    end
    return kernel
end

# The per-chunk worker for the sweep: checks every scenario in firstSignature:lastSignature for consistency, pushing the fixed points it finds onto foundFixedPoints (hence the `!` in the name).
function _sweep_chunk_all!(foundFixedPoints::Vector{Vector{Int}}, cimTranspose::Matrix{ScoreInt},
                       firstSignature::Int, lastSignature::Int,
                       variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                       numberOfDescriptors::Int, numberOfDimensions::Int,
                       margin::Int) where {ScoreInt<:Signed}
    scenario      = Vector{Int}(undef, numberOfDescriptors)  # current odometer digits
    activeRows    = Vector{Int}(undef, numberOfDescriptors)  # cimTranspose column of each descriptor's current variant
    impactBalance = zeros(ScoreInt, numberOfDimensions)

    # Decode firstSignature into its mixed-radix digits (the only division
    # in this whole chunk) and build the initial impact balance from scratch.
    remainder = firstSignature
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        variantCount = variantCounts[descriptorIndex]
        scenario[descriptorIndex] = remainder % variantCount
        activeRows[descriptorIndex] = descriptorOffsets[descriptorIndex] + scenario[descriptorIndex] + 1
        remainder = remainder ÷ variantCount
    end
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        sourceRow = activeRows[descriptorIndex]
        @simd for targetVariant in 1:numberOfDimensions
            impactBalance[targetVariant] += cimTranspose[targetVariant, sourceRow]
        end
    end

    @inbounds for currentSignature in firstSignature:lastSignature
        # Fixed-point check: is any variant strictly better (beyond the margin) than the one this scenario currently uses?
        # We only have to find one descriptor this applies to, no need to check any beyond that. So, exit at the first descriptor that fails — most scenarios fail early.
        isFixedPoint = true
        for descriptorIndex in 1:numberOfDescriptors
            offset = descriptorOffsets[descriptorIndex]

            # A challenger must beat the current variant's score PLUS the margin to unseat it.
            scoreToBeat = impactBalance[offset + scenario[descriptorIndex] + 1] + margin
            for variantColumn in 1:variantCounts[descriptorIndex]
                if impactBalance[offset + variantColumn] > scoreToBeat
                    # Some variant is strictly better, so this scenario is not consistent — no need to check the rest.
                    isFixedPoint = false
                    break #drop out of the for
                end
            end

            if !isFixedPoint
                break #drop out of the descriptors loop
            end
        end

        # Store a COPY: the scenario buffer is reused by the next iteration.
        if isFixedPoint
            push!(foundFixedPoints, copy(scenario))
        end

        # ── Advance the odometer by one and patch the impact balance ──
        # Incrementing one digit means one descriptor swapped variants, so the balance changes by (new variant's row - old variant's row): two rows, not a full recalculation.
        if currentSignature < lastSignature
            for descriptorIndex in 1:numberOfDescriptors
                variantCount = variantCounts[descriptorIndex]
                variantCount == 1 && continue    # single-variant digit never changes; carry onward
                oldRow = activeRows[descriptorIndex]
                if scenario[descriptorIndex] + 1 < variantCount
                    # Normal increment: this digit goes up by one, done.
                    scenario[descriptorIndex] += 1
                    newRow = oldRow + 1
                    activeRows[descriptorIndex] = newRow
                    @simd for targetVariant in 1:numberOfDimensions
                        impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                        cimTranspose[targetVariant, oldRow]
                    end
                    break
                end
                # This digit rolls over to 0 and the carry moves to the next descriptor (the loop continues).
                scenario[descriptorIndex] = 0
                newRow = descriptorOffsets[descriptorIndex] + 1
                activeRows[descriptorIndex] = newRow
                @simd for targetVariant in 1:numberOfDimensions
                    impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                    cimTranspose[targetVariant, oldRow]
                end
            end
        end
    end
    return foundFixedPoints
end

#endregion


#region "Branch-and-bound exhaustive search"

# ═══ How branch-and-bound works ═══════════════════════════════════════════
#
# THE PROBLEM
# -----------
# The sweep above answers "which scenarios are consistent?" by building every
# single scenario and testing it. That is honest but wasteful: a model with 10
# descriptors of 3 variants each has 59,049 scenarios, and one with 20 such
# descriptors has 3.5 billion. Usually only a handful are consistent.
#
# THE IDEA
# --------
# Instead of building whole scenarios, build them one descriptor at a time,
# and throw away whole families of them at once.
#
# Picture a tree. At the top, nothing is decided. The first level picks a
# variant for descriptor 1, the second level for descriptor 2, and so on; each
# leaf at the bottom is one complete scenario. A "node" is a partly-filled-in
# scenario, e.g. "Economy = Boom, Policy = Green, everything else undecided".
#
# The trick is that we can often prove a partial scenario is hopeless — that
# NO way of filling in the remaining descriptors could ever be consistent. When
# we can prove that, we discard the node and never visit anything beneath it.
# Cutting a node at the point where 10 descriptors of 3 variants are still
# undecided throws away 3^10 = 59,049 scenarios in a single test.
#
# HOW WE PROVE A NODE IS HOPELESS
# -------------------------------
# Remember what "consistent" means: for every descriptor, the variant the
# scenario has chosen must score at least as well as its siblings (else that
# descriptor would move, and it wouldn't be a fixed point).
#
# In a partial scenario we can't know a variant's final score, because the
# undecided descriptors will still add to it. But we can BRACKET it. A
# variant's final score is:
#
#     (what the already-decided descriptors contribute)   <- known exactly
#   + (what the undecided descriptors will contribute)    <- unknown, but bounded
#
# The first part is tracked as we descend (`prefixBalance`). For the second we
# precompute, once, the smallest and largest total the undecided descriptors
# could possibly contribute (`_bnb_bounds` → `suffixMin` / `suffixMax`).
#
# Now compare a chosen variant with one of its siblings:
#
#   * chosen variant's BEST case  = its prefix + the most the rest could add
#   * sibling's WORST case        = its prefix + the least the rest could add
#
# If the sibling's worst case still beats the chosen variant's best case, then
# the sibling wins no matter how the undecided descriptors turn out. That
# descriptor would always want to move. Every scenario below this node is
# inconsistent — prune.
#
# A WORKED EXAMPLE
# ----------------
# Three descriptors: Economy (Boom/Bust), Policy (Green/Grey), Energy
# (Renewable/Fossil). We have descended to the node
#
#     Economy = Boom, Policy = Green, Energy = still undecided
#
# and we are checking whether Policy is happy with Green. Between them, the
# descriptors decided so far contribute to Policy's two variants:
#
#     Green: -4                      Grey: +5
#
# Energy is undecided; from the precomputed bounds, it will contribute:
#
#     to Green: between  0 and +1    to Grey: between +1 and +2
#
# So whatever Energy does:
#
#     Green ends up between -4 and -3   -> its BEST possible score is  -3
#     Grey  ends up between +6 and +7   -> its WORST possible score is +6
#
# Grey's worst (+6) still beats Green's best (-3). So in every completion,
# Policy abandons Green for Grey. Both scenarios under this node — (Boom,
# Green, Renewable) and (Boom, Green, Fossil) — are inconsistent, and we
# discard them without ever scoring them. In a 12-descriptor model the same
# single test would have discarded tens of thousands.
#
# WHY THE ANSWER IS STILL EXACT
# -----------------------------
# Two things make this a shortcut rather than an approximation:
#
# 1. The bounds are conservative. We only ever prune when the sibling wins
#    even in the most favourable case for the chosen variant, so a scenario
#    that could be consistent is never discarded.
#
# 2. At the bottom of the tree, when every descriptor is decided, there is
#    nothing left undecided, so the suffix bounds are all zero and both
#    "brackets" collapse to the true score. The prune test becomes precisely
#    the ordinary consistency test. That is why a leaf which survives the
#    prune IS a fixed point and needs no further checking.
#
# The result is identical to the sweep's, right down to the ordering.
#
# WHY THERE IS A NODE BUDGET
# --------------------------
# Pruning only pays when the cross-impact matrix is strongly coupled, so that
# descriptors decisively push each other around. If the impacts are weak or
# evenly balanced, the brackets overlap almost everywhere, almost nothing gets
# pruned, and we end up walking the whole tree — which is slower than the
# sweep, because the sweep has the odometer trick and visits only the leaves
# whereas the tree walk also visits every internal node.
#
# Measured on the sample models in test/sample_files (nodes visited, as a
# percentage of the number of scenarios — lower is better):
#
#     CIB_nested         408,146,688 scenarios       33,395 nodes    0.01%
#     CIB_natl_regional   22,674,816 scenarios       27,555 nodes    0.12%
#     bench_50x50         60,466,176 scenarios    1,362,422 nodes    2.25%
#     bench_typical           59,049 scenarios        8,037 nodes   13.6%
#     bench_medium             1,024 scenarios        1,140 nodes  111%
#     CIB_global                  36 scenarios           96 nodes  267%
#
# The big, strongly-coupled models are where the win is: CIB_nested finds its
# 20 consistent scenarios after looking at one ten-thousandth of the space.
# The last two lines are the failure mode — on tiny or weakly-coupled models
# the walk costs MORE than checking every scenario, because of those internal
# nodes.
#
# We can't tell which case we're in without trying, so the search is given a
# budget of nodes it may visit. If it blows the budget, it abandons the
# attempt and `_find_kernel` falls back to the sweep. This is also why
# `:auto` doesn't even try branch-and-bound below 100,000 scenarios.
#
# ══════════════════════════════════════════════════════════════════════════

"""
    _bnb_bounds(cib) -> (suffixMin, suffixMax)

Precompute, for every variant column, how much the descriptors from `k` onwards could add to that column's score — at least (`suffixMin`) and at most (`suffixMax`). This is the "unknown but bounded" half of the bracket described in the section header, and it is computed once and reused by every node in the tree.

`suffixMin[k, c]` / `suffixMax[k, c]` is the minimum / maximum total contribution descriptors `k..numberOfDescriptors` can make to the impact score of variant column `c`, over every possible choice of their variants. Each descriptor is free to choose independently, so the extremes just add up: the minimum is the sum of each descriptor's own smallest entry in column `c`.

Row `numberOfDescriptors + 1` is left as zero — "no descriptors left to decide, so they can add nothing". That is what makes the prune test exact at the bottom of the tree. So once descriptors `1..k` are assigned, `suffixMin[k+1, c] .. suffixMax[k+1, c]` brackets what the still-undecided descriptors can add to column `c`.
"""
function _bnb_bounds(cib::CIB)
    numberOfDescriptors, numberOfDimensions = cib.numberOfDescriptors, cib.numberOfDimensions
    cimTranspose = cib.cim_t
    suffixMin = zeros(Int, numberOfDescriptors + 1, numberOfDimensions)
    suffixMax = zeros(Int, numberOfDescriptors + 1, numberOfDimensions)
    # Build the bounds back-to-front: descriptor k's bounds are its own
    # best/worst row entry plus whatever descriptors k+1..end can add.
    # `ndesc:-1:1` is a range counting down (like range(n, 0, -1) in Python).
    @inbounds for descriptorIndex in numberOfDescriptors:-1:1
        offset = cib.desc_offsets[descriptorIndex]
        variantCount = cib.numberOfVariants[descriptorIndex]
        for targetVariant in 1:numberOfDimensions
            minValue = typemax(Int)
            maxValue = typemin(Int)
            for variantIndex in 0:variantCount-1
                # = cim[row of this variant, targetVariant]
                value = Int(cimTranspose[targetVariant, offset + variantIndex + 1])
                minValue = ifelse(value < minValue, value, minValue)
                maxValue = ifelse(value > maxValue, value, maxValue)
            end
            suffixMin[descriptorIndex, targetVariant] = suffixMin[descriptorIndex + 1, targetVariant] + minValue
            suffixMax[descriptorIndex, targetVariant] = suffixMax[descriptorIndex + 1, targetVariant] + maxValue
        end
    end
    return suffixMin, suffixMax
end

# Everything one task needs while walking its part of the tree, bundled into a
# struct so the recursion passes one argument instead of a dozen.
#
# `partialScenario` is the node we are currently at — the variants decided so
# far. `prefixBalance` is what those decided descriptors contribute to every
# variant's score: the "known exactly" half of the bracket in the section
# header, kept up to date as we descend rather than recomputed at each node.
#
# The three fields after `foundFixedPoints` implement the node budget.
# `nodesVisited` is this task's private counter (a `Ref` is a mutable
# single-value box, needed because the struct itself is immutable), while
# `totalNodes` and `abortFlag` are Atomics shared by every task — the Julia
# equivalent of C#'s Interlocked operations.
struct _BnBState
    cim_t::Matrix{Int}
    sufmin::Matrix{Int}
    sufmax::Matrix{Int}
    nvariants::Vector{Int}
    offsets::Vector{Int}
    ndesc::Int
    ndim::Int
    partialScenario::Vector{Int}
    prefixBalance::Vector{Int}
    foundFixedPoints::Vector{Vector{Int}}
    nodesVisited::Base.RefValue{Int}
    totalNodes::Threads.Atomic{Int}
    abortFlag::Threads.Atomic{Bool}
    nodeBudget::Int
    margin::Int
end

"""
    _bnb_pruned(state, assignedCount) -> Bool

Can this whole branch of the tree be thrown away? This is the test worked through in the section header: with descriptors `1..assignedCount` decided, return true if some decided descriptor has a rival variant that wins in *every* possible completion — its worst case still beats the chosen variant's best case (by more than the rule's `margin`). If so, no scenario below this node can be a fixed point, so the caller skips the entire subtree.

Only *decided* descriptors are checked, because an undecided one has no chosen variant to be unhappy with yet.

Ties — and gaps within the margin — never prune, which matches the convention everywhere else that a descriptor keeps its current variant unless something strictly beats it. Once every descriptor is assigned the suffix bounds are zero, so this condition becomes exactly the margin fixed-point test on the completed scenario (`margin = 0` recovers plain global consistency).
"""
function _bnb_pruned(state::_BnBState, assignedCount::Int)
    prefixBalance = state.prefixBalance
    suffixMin = state.sufmin
    suffixMax = state.sufmax
    suffixRow = assignedCount + 1   # bounds row for the still-undecided descriptors
    @inbounds for descriptorIndex in 1:assignedCount
        offset = state.offsets[descriptorIndex]
        chosenColumn = offset + state.partialScenario[descriptorIndex] + 1
        # Best case for the variant this branch chose: what the decided
        # descriptors already give it, plus the most the undecided ones could
        # add. A rival has to clear this (and the margin) to win outright.
        # In the worked example this is Green's -4 + 1 = -3.
        bestCaseChosen = prefixBalance[chosenColumn] + suffixMax[suffixRow, chosenColumn] + state.margin
        for column in offset+1:offset+state.nvariants[descriptorIndex]
            # Worst case for a rival variant: its decided contribution plus the
            # least the undecided descriptors could add — in the example,
            # Grey's +5 + 1 = +6. If even that beats the chosen variant's best
            # case, the rival wins however the rest turns out, so nothing below
            # this node can be consistent.
            if prefixBalance[column] + suffixMin[suffixRow, column] > bestCaseChosen
                return true
            end
        end
    end
    return false
end

# Count one visited node against the budget (see "why there is a node budget"
# in the section header). Every task keeps its own private tally and only adds
# it to the shared total every 256 nodes, so the threads rarely have to touch
# shared state — the exact count doesn't matter, only whether we are roughly
# over budget. Returns false once the budget has been blown, which tells the
# whole search to give up and let the caller fall back to the sweep.
function _bnb_charge!(state::_BnBState)
    # `state.nodesVisited[]` — the [] reads/writes the value inside a Ref box.
    state.nodesVisited[] += 1
    if state.nodesVisited[] >= 256
        Threads.atomic_add!(state.totalNodes, state.nodesVisited[])
        state.nodesVisited[] = 0
        if state.totalNodes[] > state.nodeBudget || state.abortFlag[]
            state.abortFlag[] = true
            return false
        end
    end
    return true
end

# Visit one level of the tree: try each variant of descriptor `depth` in turn,
# and for each one either prune it or recurse into the descriptors below it.
# This is the depth-first walk itself — "decide one more descriptor, see if the
# result is already hopeless, and if not keep going".
#
# Only one scenario buffer exists per task, so each variant is pushed onto the
# running prefix balance before descending and popped off afterwards, leaving
# the state exactly as it was found. That is why there is no per-node
# allocation despite the recursion.
#
# Returns false if the node budget ran out and the search is abandoning.
function _bnb_node!(state::_BnBState, depth::Int)
    variantCount = state.nvariants[depth]
    offset = state.offsets[depth]
    prefixBalance = state.prefixBalance
    cimTranspose = state.cim_t
    numberOfDimensions = state.ndim
    @inbounds for variantIndex in 0:variantCount-1
        variantColumn = offset + variantIndex + 1
        # Push this variant's row onto the running prefix balance...
        @simd for targetVariant in 1:numberOfDimensions
            prefixBalance[targetVariant] += cimTranspose[targetVariant, variantColumn]
        end
        state.partialScenario[depth] = variantIndex
        keepGoing = _bnb_charge!(state)
        if keepGoing && !_bnb_pruned(state, depth)
            if depth == state.ndesc
                # All descriptors assigned and the prune (now exact) passed:
                # this is a fixed point. Store a copy — the buffer is reused.
                push!(state.foundFixedPoints, copy(state.partialScenario))
            else
                keepGoing = _bnb_node!(state, depth + 1)
            end
        end
        # ...and pop it again before trying the next variant.
        @simd for targetVariant in 1:numberOfDimensions
            prefixBalance[targetVariant] -= cimTranspose[targetVariant, variantColumn]
        end
        keepGoing || return false
    end
    return true
end

"""
    _bnb_fixed_points(cib, suffixMin, suffixMax; node_budget, margin)
        -> (Union{Nothing, Vector{Vector{Int}}}, nodesVisited)

Run the branch-and-bound search described in the section header, across all threads.

The work is split by chopping the tree near the top: every combination of variants for the first few descriptors becomes one task, which then explores its own subtree alone. Enough descriptors are taken to make comfortably more tasks than threads, so no thread sits idle waiting for a slow branch. Tasks share nothing but the node counter, and each collects its own results.

Returns a tuple. The first element is the complete kernel sorted by ascending signature — or `nothing` if the search gave up because it blew the node budget, which is the caller's signal to fall back to the sweep. The second is how many tree nodes were visited, i.e. how many partial scenarios were expanded; comparing that with the total number of scenarios shows how much the pruning actually saved.
"""
function _bnb_fixed_points(cib::CIB, sufmin::Matrix{Int}, sufmax::Matrix{Int};
                           node_budget::Int, margin::Int=0)
    numberOfDescriptors = cib.numberOfDescriptors
    variantCounts = cib.numberOfVariants
    descriptorOffsets = cib.desc_offsets
    numberOfDimensions = cib.numberOfDimensions
    cimTranspose = cib.cim_t

    # Choose how many leading descriptors to decide up front. Every combination
    # of their variants becomes one task, so keep taking descriptors until
    # there are comfortably more tasks than threads. E.g. with 4 threads and
    # 3-variant descriptors, two descriptors give 9 tasks, three give 27.
    prefixDepth = 0
    numberOfPrefixes = 1
    while prefixDepth < numberOfDescriptors && numberOfPrefixes < 4 * Threads.nthreads()
        prefixDepth += 1
        numberOfPrefixes *= variantCounts[prefixDepth]
    end

    totalNodes = Threads.Atomic{Int}(0)
    abortFlag = Threads.Atomic{Bool}(false)
    taskOutputs = [Vector{Vector{Int}}() for _ in 1:numberOfPrefixes]

    @sync for prefixId in 0:numberOfPrefixes-1
        taskOutput = taskOutputs[prefixId + 1]
        Threads.@spawn begin
            state = _BnBState(cimTranspose, sufmin, sufmax, variantCounts, descriptorOffsets,
                              numberOfDescriptors, numberOfDimensions,
                              zeros(Int, numberOfDescriptors), zeros(Int, numberOfDimensions),
                              taskOutput, Ref(0), totalNodes, abortFlag, node_budget, margin)
            # Decode this task's prefixId into variant choices for the first
            # prefixDepth descriptors (same mixed-radix arithmetic as
            # inv_signature), checking the prune at every level so a subtree
            # already dead partway through the prefix is skipped without
            # descending into it.
            remainder = prefixId
            alive = true
            @inbounds for descriptorIndex in 1:prefixDepth
                state.partialScenario[descriptorIndex] = remainder % variantCounts[descriptorIndex]
                remainder = remainder ÷ variantCounts[descriptorIndex]
                variantColumn = descriptorOffsets[descriptorIndex] +
                                state.partialScenario[descriptorIndex] + 1
                @simd for targetVariant in 1:numberOfDimensions
                    state.prefixBalance[targetVariant] += cimTranspose[targetVariant, variantColumn]
                end
                alive = _bnb_charge!(state) && !_bnb_pruned(state, descriptorIndex)
                alive || break
            end
            if alive
                if prefixDepth == numberOfDescriptors
                    # The prefix IS a full scenario and it survived the
                    # (exact, at full depth) prune: it is a fixed point.
                    push!(taskOutput, copy(state.partialScenario))
                else
                    _bnb_node!(state, prefixDepth + 1)
                end
            end
            Threads.atomic_add!(totalNodes, state.nodesVisited[])   # flush the residual private count
        end
    end

    # `abortFlag[]` reads the atomic's value.
    abortFlag[] && return (nothing, totalNodes[])
    kernel = Vector{Vector{Int}}()
    for taskOutput in taskOutputs
        append!(kernel, taskOutput)
    end
    # Tasks finish in nondeterministic order; sort so the result is stable.
    sort!(kernel; by = scenario -> signature(cib, scenario))
    return (kernel, totalNodes[])
end

#endregion

#region "basins"
"""
    find_basins(cib; rule=GlobalSuccession()) -> (fixed_points, basin_sizes, cycle_count)

Exhaustive basin-of-attraction analysis under `rule` (default [`GlobalSuccession`](@ref)). Follows the succession chain from every scenario in the space, counting how many starting points converge to each fixed point / consistent scenario.

This runs in two phases: a threaded sweep fills a flat successor table (every scenario's succession step), then a resolution pass walks the table with path compression so each scenario is resolved exactly once. Memory is ~8n bytes for the two flat tables (Int32 entries when the space fits, Int64 otherwise), independent of the thread count.
For example, the first phase works out that scenario A maps to B, B maps to C, and C maps to D which is a consistent scenario. Therefore A, B, C and D are all in D's basin, and the second phase works this out.

The default [`GlobalSuccession`](@ref) rule takes a fast path: the successor table is filled by a mixed-radix odometer that maintains the impact balance incrementally (one row delta per scenario, no allocation) instead of calling [`succession_step`](@ref) per scenario. Any other rule uses the generic per-scenario path — identical output, just slower.
"""
function find_basins(cib::CIB; rule::SuccessionRule=GlobalSuccession())
    return _basins(rule, cib)
end

# Fast path, selected by dispatch when the rule is exactly GlobalSuccession.
# The odometer chunk re-derives that rule's argmax semantics internally,
# which is exactly why this specialisation cannot serve an arbitrary rule.
function _basins(::GlobalSuccession, cib::CIB)
    numberOfScenarios = max_signature(cib) + 1
    scoreType = _score_type(cib)
    # Store table entries as Int32 when every signature fits — half the
    # memory, and memory traffic is what bounds this analysis.
    if numberOfScenarios <= Int(typemax(Int32)) - 1
        return _fast_basins(cib, Int32, Matrix{scoreType}(cib.cim_t))
    end
    return _fast_basins(cib, Int64, Matrix{scoreType}(cib.cim_t))
end

# `::Type{SignatureInt}` receives a TYPE as an argument value (Int32 or
# Int64 above) — passing types around as ordinary values is normal in Julia.
function _fast_basins(cib::CIB, ::Type{SignatureInt},
                      cimTranspose::Matrix{ScoreInt}) where {SignatureInt<:Union{Int32,Int64}, ScoreInt<:Signed}
    numberOfScenarios = max_signature(cib) + 1
    successorTable = Vector{SignatureInt}(undef, numberOfScenarios)
    _successor_table!(successorTable, cib, cimTranspose)
    fixedPointSignatures, basinSizes, cycleCount = _resolve_and_tally(successorTable, numberOfScenarios)
    fixedPoints = [inv_signature(cib, sig) for sig in fixedPointSignatures]
    return (fixedPoints, basinSizes, cycleCount)
end

# Generic path: works for any rule, because it only ever calls the rule's
# own succession_step. Allocates a few vectors per scenario, so it is much
# slower than the odometer path — a correctness baseline, not a hot loop.
function _basins(rule::SuccessionRule, cib::CIB)
    numberOfScenarios = max_signature(cib) + 1
    if numberOfScenarios <= Int(typemax(Int32)) - 1
        return _generic_basins(rule, cib, Int32)
    end
    return _generic_basins(rule, cib, Int64)
end

function _generic_basins(rule::SuccessionRule, cib::CIB,
                         ::Type{SignatureInt}) where {SignatureInt<:Union{Int32,Int64}}
    numberOfScenarios = max_signature(cib) + 1
    successorTable = Vector{SignatureInt}(undef, numberOfScenarios)
    numberOfThreads = Threads.nthreads()
    chunkSize = cld(numberOfScenarios, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * chunkSize
        lastSignature = min(threadIndex * chunkSize, numberOfScenarios) - 1
        firstSignature > lastSignature && continue
        # Each task writes a disjoint range of the table, so there is no
        # data race even though they share the array.
        Threads.@spawn for currentSignature in firstSignature:lastSignature
            scenario = inv_signature(cib, currentSignature)
            successor = succession_step(rule, cib, scenario)
            @inbounds successorTable[currentSignature + 1] =
                SignatureInt(signature(cib, successor))
        end
    end
    fixedPointSignatures, basinSizes, cycleCount = _resolve_and_tally(successorTable, numberOfScenarios)
    fixedPoints = [inv_signature(cib, sig) for sig in fixedPointSignatures]
    return (fixedPoints, basinSizes, cycleCount)
end

"""
    _successor_table!(successorTable, cib, cimTranspose) -> successorTable

Fill `successorTable[sig + 1]` with the succession-step signature of every scenario.
Threaded over contiguous chunks (disjoint writes); each chunk advances a mixed-radix odometer and maintains the impact balance incrementally, then takes the per-descriptor argmax (ties favour the current variant, then the lowest index — identical to [`succession_step`](@ref)).
"""
function _successor_table!(successorTable::Vector{SignatureInt}, cib::CIB,
                           cimTranspose::Matrix{ScoreInt}) where {SignatureInt, ScoreInt}
    numberOfScenarios = max_signature(cib) + 1
    # Mixed-radix place values, so a chunk can assemble a successor's
    # signature from its digits without any multiplication loop.
    placeValues = Vector{Int}(undef, cib.numberOfDescriptors)
    placeValue = 1
    for descriptorIndex in 1:cib.numberOfDescriptors
        placeValues[descriptorIndex] = placeValue
        placeValue *= cib.numberOfVariants[descriptorIndex]
    end

    numberOfChunks = max(1, min(numberOfScenarios, 16 * Threads.nthreads()))
    chunkSize = cld(numberOfScenarios, numberOfChunks)
    numberOfChunks = cld(numberOfScenarios, chunkSize)
    @sync for chunkIndex in 1:numberOfChunks
        firstSignature = (chunkIndex - 1) * chunkSize
        lastSignature = min(chunkIndex * chunkSize, numberOfScenarios) - 1
        Threads.@spawn _successor_chunk!(successorTable, cimTranspose,
                                         firstSignature, lastSignature,
                                         cib.numberOfVariants, cib.desc_offsets,
                                         placeValues, cib.numberOfDescriptors,
                                         cib.numberOfDimensions)
    end
    return successorTable
end

# The per-chunk worker for the successor table. Same odometer + incremental
# impact balance as _sweep_chunk_all!, but instead of a yes/no consistency
# test it records where every scenario steps to.
function _successor_chunk!(successorTable::Vector{SignatureInt}, cimTranspose::Matrix{ScoreInt},
                           firstSignature::Int, lastSignature::Int,
                           variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                           placeValues::Vector{Int}, numberOfDescriptors::Int,
                           numberOfDimensions::Int) where {SignatureInt, ScoreInt}
    scenario      = Vector{Int}(undef, numberOfDescriptors)
    activeRows    = Vector{Int}(undef, numberOfDescriptors)
    impactBalance = zeros(ScoreInt, numberOfDimensions)

    # Decode the starting signature and build the initial impact balance,
    # exactly as in _sweep_chunk_all!.
    remainder = firstSignature
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        variantCount = variantCounts[descriptorIndex]
        scenario[descriptorIndex] = remainder % variantCount
        activeRows[descriptorIndex] = descriptorOffsets[descriptorIndex] + scenario[descriptorIndex] + 1
        remainder = remainder ÷ variantCount
    end
    @inbounds for descriptorIndex in 1:numberOfDescriptors
        sourceRow = activeRows[descriptorIndex]
        @simd for targetVariant in 1:numberOfDimensions
            impactBalance[targetVariant] += cimTranspose[targetVariant, sourceRow]
        end
    end

    @inbounds for currentSignature in firstSignature:lastSignature
        # ── Successor: every descriptor independently picks its best-scoring
        #    variant, and the choices are assembled straight into a signature.
        successorSignature = 0
        for descriptorIndex in 1:numberOfDescriptors
            offset = descriptorOffsets[descriptorIndex]
            bestVariant = scenario[descriptorIndex]
            bestScore = impactBalance[offset + bestVariant + 1]  # current variant seeds the max
            for variantIndex in 0:variantCounts[descriptorIndex]-1
                score = impactBalance[offset + variantIndex + 1]
                isBetter = score > bestScore     # strict >: ties keep current / lower index
                # `ifelse` evaluates both branches and selects one — no
                # branch, so the CPU never mispredicts in this hot loop.
                bestScore = ifelse(isBetter, score, bestScore)
                bestVariant = ifelse(isBetter, variantIndex, bestVariant)
            end
            successorSignature += placeValues[descriptorIndex] * bestVariant
        end
        successorTable[currentSignature + 1] = SignatureInt(successorSignature)

        # ── Advance the odometer by one and patch the impact balance
        #    (identical to the increment in _sweep_chunk_all!).
        if currentSignature < lastSignature
            for descriptorIndex in 1:numberOfDescriptors
                variantCount = variantCounts[descriptorIndex]
                variantCount == 1 && continue    # single-variant digit never changes; carry onward
                oldRow = activeRows[descriptorIndex]
                if scenario[descriptorIndex] + 1 < variantCount
                    scenario[descriptorIndex] += 1
                    newRow = oldRow + 1
                    activeRows[descriptorIndex] = newRow
                    @simd for targetVariant in 1:numberOfDimensions
                        impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                        cimTranspose[targetVariant, oldRow]
                    end
                    break
                end
                scenario[descriptorIndex] = 0    # roll over; carry to the next digit
                newRow = descriptorOffsets[descriptorIndex] + 1
                activeRows[descriptorIndex] = newRow
                @simd for targetVariant in 1:numberOfDimensions
                    impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                    cimTranspose[targetVariant, oldRow]
                end
            end
        end
    end
    return successorTable
end

"""
    _fp_id!(registryLock, signatureForId, idForSignature, sig) -> id

Locked get-or-assign of a dense 1-based id for the fixed point with signature `sig`. Called once per fixed point discovered (≈ number-of-fixed-points times in total across all workers), so the lock is essentially uncontended. Concurrent discoverers of the same fixed point serialise here and receive the same id.
"""
# @noinline keeps this rarely-taken locked path out of the caller's hot
# loop, so the compiler optimises the loop without it.
@noinline function _fp_id!(registryLock::ReentrantLock, signatureForId::Vector{Int},
                           idForSignature::Dict{Int,Int}, fixedPointSignature::Int)
    # lock / try / finally-unlock is the standard exception-safe locking
    # pattern — the same shape as C#'s `lock` statement expands to.
    lock(registryLock)
    try
        denseId = get(idForSignature, fixedPointSignature, 0)  # 0 = not registered yet
        if denseId == 0
            push!(signatureForId, fixedPointSignature)
            denseId = length(signatureForId)
            idForSignature[fixedPointSignature] = denseId
        end
        return denseId
    finally
        unlock(registryLock)
    end
end

"""
    _resolve_chunk!(attractorLabels, successorTable, firstSignature, lastSignature,
                    registryLock, signatureForId, idForSignature)

Resolve the starting scenarios in `firstSignature:lastSignature` into the shared `attractorLabels`, following successor chains.
Cycle detection is **thread-private** (a per-worker `history` list plus a backward scan), so `attractorLabels` only ever holds *final* labels — 0 = unvisited, -1 = cycle, k > 0 = converges to the fixed point with dense id `k`. There is no shared in-progress marker, so a worker that walks into another worker's not-yet-resolved chain just re-walks it (redundant, never a false cycle) and reaches the same attractor; every scenario's attractor is deterministic, so concurrent writes to the same slot store the same value — a benign race. Fixed-point ids come from the locked registry ([`_fp_id!`](@ref)), hit only once per fixed point.
"""
function _resolve_chunk!(attractorLabels::Vector{SignatureInt}, successorTable::Vector{SignatureInt},
                         firstSignature::Int, lastSignature::Int,
                         registryLock::ReentrantLock, signatureForId::Vector{Int},
                         idForSignature::Dict{Int,Int}) where {SignatureInt}
    history = Int[]     # the chain of scenarios walked from the current start
    @inbounds for startSignature in firstSignature:lastSignature
        attractorLabels[startSignature + 1] != 0 && continue   # already resolved
        empty!(history)
        current = startSignature
        label = zero(SignatureInt)
        while true
            existingLabel = attractorLabels[current + 1]
            if existingLabel != 0
                # We walked into territory that is already resolved: the
                # whole chain behind us shares its attractor.
                label = existingLabel
                break
            end
            # Is `current` already on our own chain? Then we have walked in
            # a circle — a cycle of length ≥ 2. Scanning backwards finds a
            # repeat fastest, since a cycle closes near the chain's end.
            alreadyOnChain = false
            for position in length(history):-1:1
                if history[position] == current
                    alreadyOnChain = true
                    break
                end
            end
            if alreadyOnChain
                label = SignatureInt(-1)   # the whole chain (incl. the tail leading in) is labelled cycle
                break
            end
            push!(history, current)
            next = Int(successorTable[current + 1])
            if next == current
                # A scenario that steps to itself is a fixed point. It is
                # already in `history`, so it counts itself in its own basin.
                label = SignatureInt(_fp_id!(registryLock, signatureForId,
                                             idForSignature, current))
                break
            end
            current = next
        end
        # Path compression: every scenario on the walked chain gets the
        # final label, so later walks that touch any of them stop instantly.
        for visited in history
            attractorLabels[visited + 1] = label
        end
    end
    return nothing
end

"""
    _resolve_and_tally(successorTable, numberOfScenarios) -> (fp_sigs, sizes, cycle_count)

Resolve every scenario to its fixed point attractor by walking the successor table, then tally basin sizes and the cycle count.
"""
function _resolve_and_tally(successorTable::Vector{SignatureInt},
                            numberOfScenarios::Int) where {SignatureInt}
    attractorLabels = zeros(SignatureInt, numberOfScenarios)
    registryLock = ReentrantLock()
    signatureForId = Int[]             # dense id (1-based) -> fixed-point signature
    idForSignature = Dict{Int,Int}()   # fixed-point signature -> dense id

    numberOfThreads = Threads.nthreads()
    chunkSize = cld(numberOfScenarios, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * chunkSize
        lastSignature = min(threadIndex * chunkSize, numberOfScenarios) - 1
        Threads.@spawn _resolve_chunk!(attractorLabels, successorTable,
                                       firstSignature, lastSignature,
                                       registryLock, signatureForId, idForSignature)
    end

    # ── Threaded tally. Every label is now a dense fixed-point id (> 0) or
    #    -1 for a cycle. Each thread counts into its own arrays; the counts
    #    are combined afterwards, so no locking is needed.
    numberOfFixedPoints = length(signatureForId)
    tallyChunkSize = cld(numberOfScenarios, numberOfThreads)
    perThreadCounts = [zeros(Int, numberOfFixedPoints) for _ in 1:numberOfThreads]
    perThreadCycleCounts = zeros(Int, numberOfThreads)
    @sync for threadIndex in 1:numberOfThreads
        firstSignature = (threadIndex - 1) * tallyChunkSize
        lastSignature = min(threadIndex * tallyChunkSize, numberOfScenarios) - 1
        counts = perThreadCounts[threadIndex]
        Threads.@spawn begin
            cyclesSeen = 0
            @inbounds for signatureValue in firstSignature:lastSignature
                label = Int(attractorLabels[signatureValue + 1])
                if label == -1
                    cyclesSeen += 1
                else
                    counts[label] += 1   # label is a dense id in 1:numberOfFixedPoints
                end
            end
            perThreadCycleCounts[threadIndex] = cyclesSeen
        end
    end

    # Combine the per-thread counts into one total per fixed point.
    totalCounts = zeros(Int, numberOfFixedPoints)
    for counts in perThreadCounts
        @inbounds for fixedPointId in 1:numberOfFixedPoints
            totalCounts[fixedPointId] += counts[fixedPointId]
        end
    end
    # Discovery order depends on thread timing; sorting by signature makes
    # the output deterministic. `sortperm` returns the ordering as an index
    # list (like NumPy's argsort), applied to both arrays in step.
    sortOrder = sortperm(signatureForId)
    return signatureForId[sortOrder], totalCounts[sortOrder], sum(perThreadCycleCounts)
end

#endregion

end
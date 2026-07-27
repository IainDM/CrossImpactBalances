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

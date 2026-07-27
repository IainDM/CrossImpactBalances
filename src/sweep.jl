# ══ The sweep: check every scenario, quickly ═══════════════════════════════
#
# The obvious way to find every consistent scenario is to look at all of them.
# The obvious implementation is far too slow, though, because it repeats the
# same work over and over:
#
#   * Scoring a scenario from scratch means adding up one matrix row per
#     descriptor, for every scenario.
#   * But if the scenarios are visited in signature order, consecutive ones
#     differ in exactly ONE descriptor. Everything else about the score is
#     unchanged from the scenario before.
#
# So this sweep walks the scenarios like a car odometer — increment the last
# digit; when it runs out of variants, reset it to zero and carry into the
# next descriptor — and carries the impact balance along with it. When a digit
# ticks over, one descriptor swapped variants, so the score changes by exactly
# (the new variant's row − the old variant's row). Two rows, rather than a
# whole recomputation, and no memory allocated inside the loop at all.
#
# The work is then split across threads by giving each a contiguous stretch of
# signatures. Each thread decodes its own starting scenario once, and from
# there everything is incremental.

"""
    _find_kernel_checkall_fast(cib; margin=0) -> Vector{Vector{Int}}

Find every consistent scenario by checking all of them, using the incremental odometer walk described at the top of this file.
"""
function _find_kernel_checkall_fast(cib::CIB; margin::Int=0)
    scoreType = _score_type(cib)
    # Copy the transposed matrix into the narrowest safe integer type. The
    # inner loops are limited by how fast memory can be read, so halving the
    # width of each entry is worth the one-off copy.
    return _sweep_all_scenarios(cib, Matrix{scoreType}(cib.cim_t); margin=margin)
end

# `where {ScoreInt<:Signed}` declares a type parameter — like a C# generic
# method `Sweep<T>(...) where T : ...`. Julia compiles a separate specialised
# version of this function for each element type it is actually called with
# (Int16 or Int here). Scores can go negative, so the type must be signed.
function _sweep_all_scenarios(cib::CIB, cimTranspose::Matrix{ScoreInt};
                              margin::Int=0) where {ScoreInt<:Signed}
    numberOfScenarios = max_signature(cib) + 1

    # Split the scenarios into chunks — deliberately many more chunks than
    # threads (16 each), so that if one chunk turns out slow the other threads
    # have plenty of other work to pick up rather than sitting idle at the end.
    # `cld` is ceiling division; the min/max keep the arithmetic sane when
    # there are fewer scenarios than threads.
    numberOfChunks = max(1, min(numberOfScenarios, 16 * Threads.nthreads()))
    chunkSize = cld(numberOfScenarios, numberOfChunks)
    numberOfChunks = cld(numberOfScenarios, chunkSize)

    # Each chunk gets its own private results vector, so no two threads ever
    # write to the same place and no locking is needed anywhere.
    chunkResults = [Vector{Vector{Int}}() for _ in 1:numberOfChunks]

    # `@sync` waits for every task spawned inside the block (like
    # Task.WaitAll); `Threads.@spawn` puts one call on a worker thread (like
    # Task.Run).
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

    # The chunks covered consecutive stretches of signatures in order, so
    # simply joining them up gives a result already sorted by signature. No
    # sorting, and no duplicates to remove.
    kernel = Vector{Vector{Int}}()
    for chunkIndex in 1:numberOfChunks
        append!(kernel, chunkResults[chunkIndex])
    end
    return kernel
end

# One worker's share of the sweep: test every scenario from firstSignature to
# lastSignature, collecting the consistent ones into `found`. (The `!` in the
# name is the Julia convention warning that an argument gets modified.)
function _sweep_chunk_all!(found::Vector{Vector{Int}}, cimTranspose::Matrix{ScoreInt},
                       firstSignature::Int, lastSignature::Int,
                       variantCounts::Vector{Int}, descriptorOffsets::Vector{Int},
                       numberOfDescriptors::Int, numberOfDimensions::Int,
                       margin::Int) where {ScoreInt<:Signed}
    scenario      = Vector{Int}(undef, numberOfDescriptors)  # the odometer's current digits
    activeRows    = Vector{Int}(undef, numberOfDescriptors)  # matrix row of each descriptor's chosen variant
    impactBalance = zeros(ScoreInt, numberOfDimensions)      # carried along and patched, never rebuilt

    # Set the odometer to this chunk's starting scenario by splitting its
    # signature back into digits, then add up its impact balance from scratch.
    # This is the only division and the only full rescoring in the whole
    # chunk — everything after this point is incremental.
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
        # ── Is this scenario consistent? ──
        # Ask each descriptor in turn whether any of its variants beats the
        # one currently in use. The moment one says yes we have our answer and
        # can stop — and most scenarios fail at the first or second descriptor,
        # which is why this test is so cheap on average.
        isConsistent = true
        for descriptorIndex in 1:numberOfDescriptors
            offset = descriptorOffsets[descriptorIndex]

            # To unseat the variant in use, a rival must beat its score by
            # more than the margin (see `fixed_point_margin`). With the usual
            # margin of 0 this is just "strictly better".
            scoreToBeat = impactBalance[offset + scenario[descriptorIndex] + 1] + margin
            for variantColumn in 1:variantCounts[descriptorIndex]
                if impactBalance[offset + variantColumn] > scoreToBeat
                    # This descriptor would rather switch, so the scenario is
                    # not consistent. No point examining the others.
                    isConsistent = false
                    break
                end
            end

            if !isConsistent
                break
            end
        end

        # Store a COPY — `scenario` is the live odometer and is about to change.
        if isConsistent
            push!(found, copy(scenario))
        end

        # ── Step the odometer on to the next scenario ──
        # Adding one changes a single descriptor's variant, so the impact
        # balance shifts by (new variant's row − old variant's row) rather
        # than being recomputed. A digit that has run out of variants resets
        # to zero and the carry moves on to the next descriptor.
        if currentSignature < lastSignature
            for descriptorIndex in 1:numberOfDescriptors
                variantCount = variantCounts[descriptorIndex]
                variantCount == 1 && continue    # only one variant: nothing to change, carry onward
                oldRow = activeRows[descriptorIndex]
                if scenario[descriptorIndex] + 1 < variantCount
                    # Room left in this digit: step it up and we are done.
                    scenario[descriptorIndex] += 1
                    newRow = oldRow + 1
                    activeRows[descriptorIndex] = newRow
                    @simd for targetVariant in 1:numberOfDimensions
                        impactBalance[targetVariant] += cimTranspose[targetVariant, newRow] -
                                                        cimTranspose[targetVariant, oldRow]
                    end
                    break
                end
                # This digit is exhausted: wrap it back to its first variant
                # and let the loop carry into the next descriptor.
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
    return found
end

# Prototype: SIMD variant-major successor-table build ("vm" fast path)
#
# STATUS: parked, not wired into the package. Kept for reference.
#
# What this is
# ------------
# A drop-in replacement for the scalar per-descriptor argmax inside
# `_successor_chunk!` (src/CrossImpactBalances.jl). The scalar argmax is the
# successor-table build's hot spot (~60% of its time) and cannot auto-vectorize
# because it is a data-dependent scan over the variants of one descriptor.
# This prototype restructures the impact-balance vector into *variant-major*
# order — descriptors grouped by radix, so "variant j of 16 descriptors" is one
# contiguous 256-bit load — and argmaxes 16 descriptors at once. The blend mask
#
#     wins = (score > max) | (lane's current variant == j & score == max)
#
# reproduces the scalar strict-`>` / current-wins / lowest-index tie-break
# exactly; the successor table it produces is byte-for-byte identical.
# Measured benefit was roughly 2x on the table build on an AVX2 machine.
#
# Why it was parked
# -----------------
# It is the hardest ~160 lines in the package for a ~2x gain in one phase of
# find_basins (the resolve/tally phases are untouched by it), and it drags in
# SIMD.jl as the package's only dependency. Simplicity won.
#
# What it needs to run again
# --------------------------
# 1. `using SIMD: Vec, vload, vifelse` in the module, and SIMD.jl in
#    Project.toml ([deps] + [compat] SIMD = "3.4").
# 2. The dispatch branch in `_successor_table!` that existed at commit f533e52:
#        if S === Int32 && T === Int16 && cib.numberOfDimensions > 0
#            lay = _vm_layout(cib.numberOfVariants, cib.desc_offsets, orders, cimT)
#            ... spawn _successor_chunk_vm! per chunk ...
#        end
#    (Int16 scores are guaranteed by `_score_type`; Int32 signatures by the
#    space-size check in `find_basins`.)
# 3. The correctness harness: property_tests.jl compares find_basins output
#    against a naive chain-following oracle — run it after rewiring.
#
# Provenance: extracted verbatim from src/CrossImpactBalances.jl at f533e52.

# ─── SIMD successor chunk (Int16 scores, Int32 signatures) ──────────────────
#
# The scalar argmax above is the table build's hot spot (~60% of its time): a
# data-dependent scan of 3-4 variants per descriptor that cannot vectorize.
# This path restructures the impact-balance vector into *variant-major* order —
# grouping descriptors by radix so that "variant j of every descriptor in the
# group" is one contiguous plane — and then argmaxes 16 descriptors at once
# with 256-bit integer SIMD. Semantics are identical to the scalar path (the
# per-state successor table is byte-for-byte the same); only the evaluation
# order changes.

"""
    _vm_layout(nvariants, offsets, orders, cimT) -> NamedTuple

Build the variant-major layout for [`_successor_chunk_vm!`](@ref). Descriptors
are grouped by radix; within a group the balance slots are reordered so plane
`j` holds "variant `j` of each member descriptor" contiguously, and each group
is split into 16-lane tiles. Returns:

- `cimt_vm` — `cimT` with its slot axis permuted to variant-major (columns,
  indexed by the odometer's descriptor-major `rows`, are untouched);
- per-tile arrays `tile_r/tile_m/tile_base/tile_k0/tile_cur0` (radix, group
  size = plane stride, group base slot, lane offset, current-variant buffer
  offset);
- `ordbuf` — 16 `Int32` mixed-radix place values per tile, zero in pad lanes
  (pad lanes therefore contribute nothing to the signature);
- `curpos_of` — where each descriptor's current variant lives in the per-worker
  current-variant buffer (`ncur` entries).
"""
function _vm_layout(nvariants::Vector{Int}, offsets::Vector{Int},
                    orders::Vector{Int}, cimT::Matrix{Int16})
    ndesc = length(nvariants)
    ndim = size(cimT, 1)
    radixes = sort!(unique(nvariants))
    dm_of_vm = Vector{Int}(undef, ndim)      # variant-major slot -> descriptor-major slot
    curpos_of = zeros(Int, ndesc)
    tile_r = Int[]; tile_m = Int[]; tile_base = Int[]; tile_k0 = Int[]; tile_cur0 = Int[]
    ordbuf = Int32[]
    base = 0
    for r in radixes
        members = [i for i in 1:ndesc if nvariants[i] == r]
        m = length(members)
        for j in 0:r-1, (k, i) in enumerate(members)
            dm_of_vm[base + j*m + k] = offsets[i] + j + 1
        end
        for k0 in 0:16:m-1
            lanes = min(16, m - k0)
            push!(tile_r, r); push!(tile_m, m); push!(tile_base, base); push!(tile_k0, k0)
            push!(tile_cur0, length(tile_r) * 16 - 15)     # 16 cur lanes per tile
            for l in 1:16
                push!(ordbuf, l <= lanes ? Int32(orders[members[k0 + l]]) : Int32(0))
            end
            for l in 1:lanes
                curpos_of[members[k0 + l]] = (length(tile_r) - 1) * 16 + l
            end
        end
        base += r * m
    end
    return (cimt_vm = cimT[dm_of_vm, :], ntiles = length(tile_r),
            tile_r = tile_r, tile_m = tile_m, tile_base = tile_base,
            tile_k0 = tile_k0, tile_cur0 = tile_cur0, ordbuf = ordbuf,
            curpos_of = curpos_of, ncur = length(tile_r) * 16)
end

"""
    _successor_chunk_vm!(succ, lay, first_sig, last_sig, nvariants, offsets,
                         ndesc, ndim)

Variant-major successor chunk: identical output to [`_successor_chunk!`](@ref)
(same incremental odometer, same strict-`>`/current-wins tie-break), with the
per-descriptor argmax vectorized across each 16-lane tile. Per variant plane
the winner rule is a single blend:

    wins = (score > max) | (lane's current variant == j  &  score == max)

which reproduces the scalar tie-break exactly: the current variant beats an
equal incumbent (the `==` arm fires only on the lane's own current plane), all
other variants need strict `>`, and lower indices win otherwise because planes
are scanned in ascending order. Pad lanes are harmless: their loads stay within
the 16-slot padding of `ibv` and their place values in `ordbuf` are zero. The
signature is then one SIMD widening multiply + horizontal sum per tile.
"""
function _successor_chunk_vm!(succ::Vector{Int32}, lay, first_sig::Int, last_sig::Int,
                              nvariants::Vector{Int}, offsets::Vector{Int},
                              ndesc::Int, ndim::Int)
    cimt_vm = lay.cimt_vm
    v    = Vector{Int}(undef, ndesc)
    rows = Vector{Int}(undef, ndesc)
    ibv  = zeros(Int16, ndim + 16)           # +16: tile loads may overhang the last group
    cur  = zeros(Int16, lay.ncur)

    s = first_sig
    @inbounds for i in 1:ndesc
        nv = nvariants[i]
        v[i] = s % nv
        rows[i] = offsets[i] + v[i] + 1
        cur[lay.curpos_of[i]] = Int16(v[i])
        s = s ÷ nv
    end
    @inbounds for i in 1:ndesc
        r = rows[i]
        @simd for j in 1:ndim
            ibv[j] += cimt_vm[j, r]
        end
    end

    V = Vec{16,Int16}
    W = Vec{16,Int32}
    @inbounds for sig in first_sig:last_sig
        # ── Successor signature: tiled SIMD argmax over the variant planes ──
        w_sig = 0
        for t in 1:lay.ntiles
            r = lay.tile_r[t]; m = lay.tile_m[t]
            off0 = lay.tile_base[t] + lay.tile_k0[t] + 1
            cv = vload(V, cur, lay.tile_cur0[t])
            mx = vload(V, ibv, off0)                     # plane 0 seeds (lowest index)
            best = V(Int16(0))
            for j in 1:r-1
                pj = vload(V, ibv, off0 + j*m)
                jv = V(Int16(j))
                b = (pj > mx) | ((cv == jv) & (pj == mx))
                mx = vifelse(b, pj, mx)
                best = vifelse(b, jv, best)
            end
            w_sig += Int(sum(convert(W, best) * vload(W, lay.ordbuf, (t-1)*16 + 1)))
        end
        succ[sig + 1] = Int32(w_sig)

        # ── Odometer increment with fused row-delta ibv update ──
        if sig < last_sig
            for i in 1:ndesc
                nv = nvariants[i]
                nv == 1 && continue      # radix-1: value stays 0, carry onward
                rold = rows[i]
                if v[i] + 1 < nv
                    vi = v[i] + 1
                    v[i] = vi
                    rnew = rold + 1
                    rows[i] = rnew
                    cur[lay.curpos_of[i]] = Int16(vi)
                    @simd for j in 1:ndim
                        ibv[j] += cimt_vm[j, rnew] - cimt_vm[j, rold]
                    end
                    break
                end
                v[i] = 0                 # roll over; carry to next digit
                rnew = offsets[i] + 1
                rows[i] = rnew
                cur[lay.curpos_of[i]] = Int16(0)
                @simd for j in 1:ndim
                    ibv[j] += cimt_vm[j, rnew] - cimt_vm[j, rold]
                end
            end
        end
    end
    return succ
end

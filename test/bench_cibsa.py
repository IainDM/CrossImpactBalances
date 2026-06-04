#!/usr/bin/env python3
"""Three-way CIB benchmark — CIBSA (upstream sei-international/cibsa) leg.

Imports the cloned package at D:\\GitHub\\cibsa (minimally ported to Python 3 —
see that repo's `git diff`). Times find_consistent on a PRE-PARSED model
(kernel=[] skips the constructor search), so parse cost is excluded and the
numbers compare like-for-like with JuCIB find-only.

Two modes per file:
  - full enum  (mc_threshold=10**18)  — apples-to-apples with JuCIB full-enum-1T.
  - MC sample  (mc_threshold=10000, seeded) — for spaces too large to enumerate.

Run:  python test/bench_cibsa.py
"""
import sys, os, json, time

CIBSA_DIR = r"D:\GitHub\cibsa"
sys.path.insert(0, CIBSA_DIR)
import numpy as np
from CIB_sim_anneal import CIB

HERE   = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, "sample_files")
MC_THRESHOLD = 10000
MC_SEED = 999
FULL = 10 ** 18           # force full enumeration

# (name, run_full_enum, run_mc)  — full enum is infeasible in pure Python on 60M.
SPECS = [
    ("CIB_global",    True,  False),
    ("bench_medium",  True,  False),
    ("bench_large",   True,  False),
    ("bench_typical", True,  True),
    ("bench_xlarge",  True,  True),
    ("bench_50x50",   False, True),
]


def median(ts):
    ts = sorted(ts)
    return ts[len(ts) // 2] if len(ts) % 2 else ts[len(ts) // 2 - 1]


def parse_seconds(path, n):
    CIB(path, kernel=[], mc_threshold=FULL)          # warmup
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        CIB(path, kernel=[], mc_threshold=FULL)
        ts.append(time.perf_counter() - t0)
    return median(ts)


def find_seconds(path, mc, n, seed=None):
    """Time find_consistent on a pre-parsed model; return (median_s, sorted_sigs, size)."""
    base = CIB(path, kernel=[], mc_threshold=mc)     # parse only
    if seed is not None:
        np.random.seed(seed)
    base.find_consistent()                           # warmup
    ts, kern = [], []
    for _ in range(n):
        if seed is not None:
            np.random.seed(seed)
        t0 = time.perf_counter()
        kern = base.find_consistent()
        ts.append(time.perf_counter() - t0)
    sigs = sorted(int(base.signature(u)) for u in kern)
    return median(ts), sigs, len(kern)


def main():
    print(f"CIBSA benchmark — Python {sys.version.split()[0]}, numpy {np.__version__}")
    print(f"package: {CIBSA_DIR}")
    results = []
    for name, do_full, do_mc in SPECS:
        path = os.path.join(SAMPLE, f"{name}.scw")
        if not os.path.isfile(path):
            print(f"  MISSING: {path} — skipping")
            continue
        probe = CIB(path, kernel=[], mc_threshold=FULL)
        total = probe.max_signature() + 1
        big = total > 1_000_000
        n = 2 if big else 3
        print(f"\n=== {name}  ({total} scenarios, {probe.ndesc} descriptors) ===")

        parse_s = parse_seconds(path, n)

        full_s = full_sigs = None
        full_size = None
        if do_full:
            full_s, full_sigs, full_size = find_seconds(path, FULL, n)
            print(f"  full-enum: {full_s:.4f}s  kernel={full_size}  sigs={full_sigs}")
        else:
            print("  full-enum: DNF (skipped — pure-Python enumeration of "
                  f"{total} scenarios is infeasible)")

        mc_s = mc_sigs = None
        mc_size = None
        if do_mc:
            mc_s, mc_sigs, mc_size = find_seconds(path, MC_THRESHOLD, n, seed=MC_SEED)
            print(f"  MC({MC_THRESHOLD}): {mc_s:.4f}s  kernel={mc_size}  sigs={mc_sigs}")

        results.append({
            "file": name, "tool": "cibsa", "total_scenarios": int(total),
            "parse_s": round(parse_s, 6),
            "full_enum_s": None if full_s is None else round(full_s, 6),
            "full_enum_kernel_size": full_size,
            "full_enum_kernel_sigs": full_sigs,
            "mc_s": None if mc_s is None else round(mc_s, 6),
            "mc_threshold": MC_THRESHOLD,
            "mc_seed": MC_SEED,
            "mc_kernel_size": mc_size,
            "mc_kernel_sigs": mc_sigs,
            "notes": "full-enum DNF (60M, pure Python)" if not do_full else "",
        })

    outpath = os.path.join(HERE, "bench_results_cibsa.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2,
                  default=lambda o: int(o) if hasattr(o, "__int__") else str(o))
    print(f"\nWrote {outpath}  ({len(results)} files)")


if __name__ == "__main__":
    main()

"""
Generate a 50x50 CIM (15 descriptors) and run Python CIBSA.
Save the .scw file and results for Julia comparison.

15 descriptors: 5 with 4 variants + 10 with 3 variants = 50 CIM dimension
Scenario space: 4^5 * 3^10 = 60,466,176 (deep MC sampling)
"""
import numpy as np
import json
import time
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from generate_benchmarks import CIB, generate_scw, NumpyEncoder

SAMPLE_DIR = os.path.join(os.path.dirname(__file__), "sample_files")

nvars = [4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
ndesc = 15
ndim = sum(nvars)
total = int(np.prod(nvars))

print(f"50x50 CIM: {ndesc} descriptors, variants={nvars}")
print(f"CIM dimension: {ndim}x{ndim}, scenario space: {total:,}")

# Use seed=0 for generating the CIM, and fix numpy state before each
# CIB run so MC sampling is reproducible
scw_path = os.path.join(SAMPLE_DIR, "bench_50x50.scw")
generate_scw(scw_path, ndesc, nvars, value_range=2, seed=0)

# Benchmark: find_consistent (3 runs, take median)
# Fix random state before each run for reproducible MC sampling
print(f"\nRunning Python CIBSA (3 trials)...")
times_find = []
for trial in range(3):
    np.random.seed(999)  # fix MC sampling seed
    t0 = time.perf_counter()
    cib = CIB(scw_path)
    elapsed = time.perf_counter() - t0
    times_find.append(elapsed)
    print(f"  Trial {trial+1}: {elapsed:.4f}s, {len(cib.kernel)} consistent")

t_find = sorted(times_find)[1]  # median

# Sort kernel and verify
cib.kernel.sort(key=lambda u: cib.signature(u))
kernel_sigs = [cib.signature(u) for u in cib.kernel]

print(f"\nVerifying {len(cib.kernel)} consistent scenarios are true fixed points...")
for u in cib.kernel:
    v = cib.succession_step(u)
    assert u == v, f"  FAIL: {u} -> {v}"
print("  All verified.")

# Benchmark IPM
times_ipm = []
for trial in range(3):
    t0 = time.perf_counter()
    ipm = cib.inner_product_matrix()
    times_ipm.append(time.perf_counter() - t0)
t_ipm = sorted(times_ipm)[1]

print(f"\nResults:")
print(f"  Kernel size: {len(cib.kernel)}")
print(f"  Kernel sigs: {kernel_sigs}")
print(f"  find_consistent: {t_find:.4f}s (median)")
print(f"  inner_product_matrix ({len(cib.kernel)}x{len(cib.kernel)}): {t_ipm:.6f}s")

# Save — include the actual kernel vectors so Julia can verify them as fixed points
results = {
    "nvariants": nvars,
    "total_scenarios": total,
    "seed": 0,
    "value_range": 2,
    "kernel_sigs": kernel_sigs,
    "kernel": [list(u) for u in cib.kernel],
    "inner_product_matrix": ipm,
    "n_kernel": len(cib.kernel),
    "python_time_find_consistent_s": round(t_find, 6),
    "python_time_ipm_s": round(t_ipm, 6),
}

outpath = os.path.join(SAMPLE_DIR, "bench_50x50_expected.json")
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2, cls=NumpyEncoder)

print(f"\nSaved to {outpath}")

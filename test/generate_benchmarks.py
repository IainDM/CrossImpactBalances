"""
Generate random .scw files and run Python CIBSA on them.
Outputs JSON with kernel signatures, inner product matrices, and timings.
"""
import numpy as np
import re
import json
import time
import os

# ── Reusable CIB class (Python CIBSA port) ──────────────────────────────────

class CIB:
    def __init__(self, scw_file, sl_file=None, kernel=None, mc_threshold=10000):
        self.descriptors = []
        self.variants = {}
        self.nvariants = []
        self.mc_threshold = mc_threshold
        s = 0; n = 0; r = 0; d = -1; v = 0
        current_desc = ""
        cim_rows = []
        with open(scw_file, 'r') as f:
            for line in f:
                stripped = line.lstrip()
                if not stripped or stripped == '\n':
                    continue
                if s == 0:
                    if stripped[0] == '&':
                        desc = stripped[1:].strip()
                        self.descriptors.append(desc)
                        self.variants[desc] = []
                        if d > -1:
                            self.nvariants.append(v)
                        v = 0; d += 1; current_desc = desc
                    elif stripped[0] == '-':
                        self.variants[current_desc].append(stripped[1:].strip())
                        n += 1; v += 1
                    elif stripped[0] == '#':
                        s = 1; self.nvariants.append(v)
                elif s < 5 and stripped[0] == '#':
                    s += 1
                elif s == 5:
                    if stripped[0] == '#':
                        s += 1
                    else:
                        cim_rows.append(list(map(int, stripped.split(','))))
                        r += 1

        self.cim = np.array(cim_rows)
        self.ndim = n
        self.ndesc = d + 1
        self.thresholds = [0] * self.ndesc

        if kernel is not None:
            self.kernel = kernel
        elif sl_file is not None:
            self.kernel = self._load_sl(sl_file)
        else:
            self.kernel = self.find_consistent()

    def _load_sl(self, sl_file):
        kernel = []
        pat = re.compile(r'^"([^"]+)"')
        with open(sl_file, 'r') as f:
            for line in f:
                stripped = line.lstrip()
                if not stripped or stripped[0] != '"':
                    continue
                m = pat.match(stripped)
                if m:
                    kernel.append([x - 1 for x in map(int, m.group(1).split())])
        return kernel

    def _varndx_to_tablendx(self, u):
        ndx = list(u)
        offset = 0
        for i in range(len(u)):
            ndx[i] = offset + u[i]
            offset += self.nvariants[i]
        return ndx

    def impact_balance(self, u):
        return np.sum(self.cim[self._varndx_to_tablendx(u)], axis=0)

    def own_impact_balance(self, u):
        return self.impact_balance(u)[self._varndx_to_tablendx(u)]

    def cross_impact_balance(self, u, v):
        return self.impact_balance(u)[self._varndx_to_tablendx(v)]

    def inner_product(self, u, v):
        ib = self.impact_balance(u)
        return int(sum(ib[self._varndx_to_tablendx(v)]))

    def succession_step(self, u):
        ib = self.impact_balance(u)
        v = list(u)
        start = 0
        for i in range(self.ndesc):
            stop = start + self.nvariants[i]
            ib_desc = ib[start:stop]
            max_val = ib_desc[u[i]]
            for j in range(len(ib_desc)):
                if ib_desc[j] > max_val:
                    max_val = ib_desc[j]
                    v[i] = j
            start = stop
        return v

    def signature(self, u):
        order = 1; sig = 0
        for ui, nv in zip(u, self.nvariants):
            sig += order * ui
            order *= nv
        return sig

    def inv_signature(self, s):
        u = []
        for nv in self.nvariants:
            u.append(s % nv)
            s //= nv
        return u

    def max_signature(self):
        return self.signature([n - 1 for n in self.nvariants])

    def get_scenario_signatures(self):
        n = self.max_signature() + 1
        if n > self.mc_threshold:
            return np.random.choice(range(n), self.mc_threshold, replace=False)
        return range(n)

    def succession(self, u):
        iterations_sig = [self.signature(u)]
        v = list(u)
        while True:
            v = self.succession_step(v)
            v_sig = self.signature(v)
            n = 1; foundit = False
            for hist in iterations_sig[::-1]:
                if hist == v_sig:
                    foundit = True; break
                n += 1
            if foundit:
                return [n, v]
            iterations_sig.append(v_sig)

    def find_consistent(self, ignore_cycles=True):
        kernel = []; signatures = set()
        for v_sig in self.get_scenario_signatures():
            v = self.inv_signature(v_sig)
            nper, veqm = self.succession(v)
            if ignore_cycles and nper > 1:
                continue
            veqm_sig = self.signature(veqm)
            if veqm_sig not in signatures:
                signatures.add(veqm_sig)
                kernel.append(veqm)
        return kernel

    def inner_product_matrix(self):
        M = []
        for u in self.kernel:
            row = []
            for v in self.kernel:
                row.append(self.inner_product(u, v))
            M.append(row)
        return M


# ── Generate random .scw file ───────────────────────────────────────────────

def generate_scw(filename, ndesc, nvars, value_range=6, seed=None):
    """Generate a random .scw file with given structure."""
    if seed is not None:
        np.random.seed(seed)

    ndim = sum(nvars)

    # Generate random CIM with zeros on diagonal blocks
    cim = np.random.randint(-value_range, value_range + 1, size=(ndim, ndim))

    # Zero out diagonal blocks
    offset = 0
    for nv in nvars:
        cim[offset:offset+nv, offset:offset+nv] = 0
        offset += nv

    with open(filename, 'w') as f:
        f.write("$ ScenarioWizard 4.0\n")
        f.write(f"{os.path.basename(filename)}\n")
        for i in range(ndesc):
            f.write(f"&D{i+1}\n")
            for j in range(nvars[i]):
                f.write(f" -V{i+1}_{j+1}\n")
        f.write("#\n")
        for i in range(ndesc):
            f.write("&\n")
            for j in range(nvars[i]):
                f.write(" \n")
        f.write("#\n#\n")
        for i in range(ndim):
            f.write("FFFFFF\n")
        f.write("#\n#\n")
        for i in range(ndim):
            row_str = ",".join(f"{v:2d} " for v in cim[i])
            f.write(row_str + "\n")
        f.write("#\n# 0-0\n\"\"\n#\n#\n")


# ── Scan seeds to find ones that produce consistent scenarios ────────────────

def find_good_seed(ndesc, nvars, value_range=3, max_seed=200):
    """Find a seed that produces at least 2 consistent scenarios."""
    for seed in range(max_seed):
        np.random.seed(seed)
        ndim = sum(nvars)
        cim = np.random.randint(-value_range, value_range + 1, size=(ndim, ndim))
        offset = 0
        for nv in nvars:
            cim[offset:offset+nv, offset:offset+nv] = 0
            offset += nv
        # Quick check: try a few starting points
        # Write temp file
        tmpf = "/tmp/_seed_test.scw"
        generate_scw(tmpf, ndesc, nvars, value_range, seed)
        try:
            c = CIB(tmpf)
            if len(c.kernel) >= 2:
                return seed, len(c.kernel)
        except:
            pass
    return None, 0


# ── Main ─────────────────────────────────────────────────────────────────────

SAMPLE_DIR = "/home/user/CrossImpactBalances/test/sample_files"
results = {}

# Medium: 5 descriptors, 4 variants each = 1024 scenarios
# Use lower value_range to increase chance of fixed points
print("Scanning for good seeds (medium)...")
seed_med, nk_med = find_good_seed(5, [4, 4, 4, 4, 4], value_range=3)
print(f"  Medium: seed={seed_med}, {nk_med} consistent scenarios")

# Large: 8 descriptors, 3 variants each = 6561 scenarios
print("Scanning for good seeds (large)...")
seed_lrg, nk_lrg = find_good_seed(8, [3, 3, 3, 3, 3, 3, 3, 3], value_range=3)
print(f"  Large: seed={seed_lrg}, {nk_lrg} consistent scenarios")

# Extra-large: 10 descriptors, 3 variants each = 59049 scenarios (uses MC sampling)
print("Scanning for good seeds (xlarge, MC sampling)...")
seed_xl, nk_xl = find_good_seed(10, [3, 3, 3, 3, 3, 3, 3, 3, 3, 3], value_range=2, max_seed=50)
print(f"  XLarge: seed={seed_xl}, {nk_xl} consistent scenarios")

configs = [
    ("bench_medium", 5, [4, 4, 4, 4, 4], 3, seed_med),
    ("bench_large",  8, [3, 3, 3, 3, 3, 3, 3, 3], 3, seed_lrg),
]
if seed_xl is not None:
    configs.append(("bench_xlarge", 10, [3, 3, 3, 3, 3, 3, 3, 3, 3, 3], 2, seed_xl))

for name, ndesc, nvars, vrange, seed in configs:
    scw_path = os.path.join(SAMPLE_DIR, f"{name}.scw")
    print(f"\n{'='*60}")
    print(f"{name}: {ndesc} descriptors, variants={nvars}, "
          f"total={np.prod(nvars)}, seed={seed}")

    generate_scw(scw_path, ndesc, nvars, vrange, seed)

    # Benchmark: find_consistent (3 runs, take median)
    times_find = []
    for trial in range(3):
        t0 = time.perf_counter()
        cib = CIB(scw_path)
        times_find.append(time.perf_counter() - t0)

    # Sort kernel by signature for deterministic ordering
    cib.kernel.sort(key=lambda u: cib.signature(u))
    kernel_sigs = [cib.signature(u) for u in cib.kernel]
    t_find = sorted(times_find)[1]  # median of 3
    print(f"  {len(cib.kernel)} consistent scenarios")
    print(f"  find_consistent: {t_find:.4f}s (median of 3: {times_find})")

    # Benchmark: inner_product_matrix (3 runs)
    times_ipm = []
    for trial in range(3):
        t0 = time.perf_counter()
        ipm = cib.inner_product_matrix()
        times_ipm.append(time.perf_counter() - t0)

    t_ipm = sorted(times_ipm)[1]
    print(f"  inner_product_matrix: {t_ipm:.6f}s (median of 3)")

    results[name] = {
        "nvariants": nvars,
        "total_scenarios": int(np.prod(nvars)),
        "seed": seed,
        "value_range": vrange,
        "kernel_sigs": kernel_sigs,
        "inner_product_matrix": ipm,
        "n_kernel": len(cib.kernel),
        "python_time_find_consistent_s": round(t_find, 6),
        "python_time_ipm_s": round(t_ipm, 6),
    }

# Save (convert numpy ints to native Python ints)
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

outpath = os.path.join(SAMPLE_DIR, "benchmark_expected.json")
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2, cls=NumpyEncoder)

print(f"\n{'='*60}")
print(f"Results saved to {outpath}")
for name, r in results.items():
    print(f"  {name}: {r['total_scenarios']} scenarios, "
          f"{r['n_kernel']} consistent, "
          f"python find={r['python_time_find_consistent_s']:.4f}s, "
          f"python ipm={r['python_time_ipm_s']:.6f}s")

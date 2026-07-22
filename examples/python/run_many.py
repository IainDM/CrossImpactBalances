"""Python interface — batch analysis and in-place editing.

Runs the same engine as the Julia examples, driven from Python via juliacall.
Install the package first (from the repo root):

    pip install -e "python/[pandas]"

Then, for multi-threaded exhaustive/basin analysis:

    PYTHON_JULIACALL_THREADS=auto python examples/python/run_many.py
"""

import os

from crossimpactbalances import Model, run_models, sweep_impact

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DIR = os.path.join(HERE, "..", "..", "test", "sample_files")


def main():
    # 1. Batch over every sample .scw file: basin analysis, one row per model.
    print("=== Batch basin analysis over sample models ===")
    records = run_models(SAMPLE_DIR, analysis="basins")
    for rec in records:
        name = os.path.basename(rec["path"])
        if "error" in rec:
            print(f"  {name:40s} ERROR: {rec['error']}")
            continue
        print(f"  {name:40s} "
              f"{rec['n_descriptors']:>3} descriptors, "
              f"{rec['n_scenarios']:>12,} scenarios, "
              f"{rec['n_consistent']:>3} consistent, "
              f"largest basin {rec['largest_basin_fraction']:.1%}")

    # 2. Load one model once, then edit values in place and re-run — no re-parse.
    print("\n=== In-place editing on CIB_global ===")
    m = Model.load(os.path.join(SAMPLE_DIR, "CIB_global.scw"))
    print(f"  {m}")
    print(f"  consistent scenarios (baseline): {len(m.consistent_scenarios())}")

    d1, d2 = m.descriptors[0], m.descriptors[1]
    src = (d1, m.variants[d1][0])
    tgt = (d2, m.variants[d2][0])
    print(f"  sweeping impact {src} -> {tgt} without reloading the file:")
    for rec in sweep_impact(m, src, tgt, [-6, -3, 0, 3, 6], analysis="consistent"):
        print(f"    value={rec['value']:>3}  ->  {rec['n_consistent']} consistent scenarios")


if __name__ == "__main__":
    main()

"""Integration tests for the batch helpers (require the Julia engine)."""

import os

import pytest

pytestmark = pytest.mark.usefixtures("engine_available")


def test_run_models_over_directory(sample_dir):
    from crossimpactbalances import run_models
    # A few small, fast sample files.
    paths = [os.path.join(sample_dir, f)
             for f in ("CIB_global.scw", "CIB_natl_regional.scw")]
    records = run_models(paths, analysis="consistent")
    assert len(records) == len(paths)
    for rec in records:
        assert "error" not in rec, rec
        assert rec["n_consistent"] == len(rec["scenarios"])
        assert rec["n_descriptors"] >= 1


def test_run_models_basins_invariant(sample_dir):
    from crossimpactbalances import run_models
    path = os.path.join(sample_dir, "CIB_global.scw")
    (rec,) = run_models([path], analysis="basins")
    covered = sum(r["basin_size"] for r in rec["scenarios"]) + rec["cycle_count"]
    assert covered == rec["n_scenarios"]


def test_run_models_records_errors(tmp_path):
    from crossimpactbalances import run_models
    bad = tmp_path / "missing.scw"
    records = run_models([str(bad)], analysis="consistent")
    assert len(records) == 1
    assert "error" in records[0]


def test_sweep_impact_restores(model):
    from crossimpactbalances import sweep_impact
    m = model
    d1, d2 = m.descriptors[0], m.descriptors[1]
    src, tgt = (d1, m.variants[d1][0]), (d2, m.variants[d2][0])
    original = m.get_impact(src, tgt)

    records = sweep_impact(m, src, tgt, [-3, 0, 3, 6], analysis="consistent")
    assert [r["value"] for r in records] == [-3, 0, 3, 6]
    assert all("n_consistent" in r for r in records)
    # Model restored to its original value afterwards.
    assert m.get_impact(src, tgt) == original

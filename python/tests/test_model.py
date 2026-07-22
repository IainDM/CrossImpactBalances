"""Integration tests for the Model wrapper (require the Julia engine)."""

import os

import pytest

pytestmark = pytest.mark.usefixtures("engine_available")


def test_structure(global_model):
    m = global_model
    assert m.n_descriptors == 3
    assert set(m.variants) == set(m.descriptors)
    # CIB_global is a 3-descriptor, 36-scenario space.
    assert m.n_scenarios == 36


def test_consistent_scenarios_are_name_maps(global_model):
    ks = global_model.consistent_scenarios()
    assert ks, "expected at least one consistent scenario"
    for sc in ks:
        assert isinstance(sc, dict)
        assert set(sc) == set(global_model.descriptors)
        for desc, variant in sc.items():
            assert variant in global_model.variants[desc]


def test_consistent_matches_loaded_kernel(global_model):
    """The computed kernel signatures match the ones in the .sl file."""
    from crossimpactbalances._engine import get_engine
    eng = get_engine()
    loaded = {int(eng.signature(global_model._cib, u))
              for u in global_model._cib.kernel}
    computed = {global_model.signature(sc)
                for sc in global_model.consistent_scenarios()}
    assert computed == loaded


def test_signature_roundtrip(global_model):
    for sc in global_model.consistent_scenarios():
        s = global_model.signature(sc)
        assert global_model.scenario_from_signature(s) == sc


def test_impact_balance_and_consistency(global_model):
    sc = global_model.consistent_scenarios()[0]
    ib = global_model.impact_balance(sc)
    assert ib["is_consistent"] is True
    assert global_model.is_consistent(sc) is True
    assert set(ib["scores"]) == set(global_model.descriptors)


def test_succession_converges_from_consistent(global_model):
    sc = global_model.consistent_scenarios()[0]
    trace = global_model.succession(sc)
    assert trace["converged"] is True
    assert trace["cycle_length"] == 1
    assert trace["steps"][0] == sc


def test_basins_cover_space(global_model):
    result = global_model.basins()
    assert result.total == global_model.n_scenarios
    covered = sum(r["basin_size"] for r in result.scenarios) + result.cycle_count
    assert covered == result.total
    # Ranked descending.
    sizes = [r["basin_size"] for r in result.scenarios]
    assert sizes == sorted(sizes, reverse=True)


def test_set_impact_edits_in_place(global_model):
    m = global_model
    d1, d2 = m.descriptors[0], m.descriptors[1]
    v1 = m.variants[d1][0]
    v2 = m.variants[d2][1]
    src, tgt = (d1, v1), (d2, v2)

    original = m.get_impact(src, tgt)
    prev = m.set_impact(src, tgt, original + 7)
    assert prev == original
    assert m.get_impact(src, tgt) == original + 7
    # Index form resolves to the same cell.
    assert m.get_impact((0, 0), (1, 1)) == original + 7
    m.set_impact(src, tgt, original)
    assert m.get_impact(src, tgt) == original


def test_edit_changes_consistent_set_without_reload(global_model):
    m = global_model
    d1, d2 = m.descriptors[0], m.descriptors[1]
    before = {m.signature(sc) for sc in m.consistent_scenarios()}

    # Crank a couple of cells hard to move the fixed points.
    m.set_impact((d2, m.variants[d2][0]), (d1, m.variants[d1][0]), 50)
    m.set_impact((d2, m.variants[d2][-1]), (d1, m.variants[d1][0]), -50)
    after = {m.signature(sc) for sc in m.consistent_scenarios()}
    assert after != before


def test_copy_is_independent(global_model):
    m = global_model
    d1, d2 = m.descriptors[0], m.descriptors[1]
    src, tgt = (d1, m.variants[d1][0]), (d2, m.variants[d2][0])
    baseline = m.get_impact(src, tgt)

    branch = m.copy()
    branch.set_impact(src, tgt, baseline + 100)
    assert branch.get_impact(src, tgt) == baseline + 100
    assert m.get_impact(src, tgt) == baseline   # original untouched


def test_bad_scenario_raises(global_model):
    with pytest.raises(ValueError):
        global_model.impact_balance({"nope": "x"})

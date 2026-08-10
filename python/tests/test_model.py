"""Integration tests for the Model wrapper. Backend-agnostic (run once per
backend via CIB_TEST_BACKEND)."""

import pytest

pytestmark = pytest.mark.usefixtures("engine_available")


def test_structure(model):
    assert model.n_descriptors == 3
    assert set(model.variants) == set(model.descriptors)
    assert model.n_scenarios == 36            # 3×3×4 space


def test_consistent_scenarios_are_name_maps(model):
    ks = model.consistent_scenarios()
    assert ks, "expected at least one consistent scenario"
    for sc in ks:
        assert isinstance(sc, dict)
        assert set(sc) == set(model.descriptors)
        for desc, variant in sc.items():
            assert variant in model.variants[desc]


def test_consistent_are_self_consistent(model):
    """Every scenario returned as consistent must be a fixed point."""
    ks = model.consistent_scenarios()
    for sc in ks:
        assert model.is_consistent(sc) is True


def test_exhaustive_matches_default(model):
    a = {model.signature(sc) for sc in model.consistent_scenarios()}
    b = {model.signature(sc) for sc in model.consistent_scenarios(exhaustive=True)}
    assert a == b


def test_signature_roundtrip(model):
    for sc in model.consistent_scenarios():
        s = model.signature(sc)
        assert model.scenario_from_signature(s) == sc


def test_impact_balance_and_consistency(model):
    sc = model.consistent_scenarios()[0]
    ib = model.impact_balance(sc)
    assert ib["is_consistent"] is True
    assert set(ib["scores"]) == set(model.descriptors)


def test_succession_converges_from_consistent(model):
    sc = model.consistent_scenarios()[0]
    trace = model.succession(sc)
    assert trace["converged"] is True
    assert trace["cycle_length"] == 1
    assert trace["steps"][0] == sc


def test_succession_reports_a_real_cycle(model):
    """A trajectory that ends in a cycle must not read as converged.

    All-first-variants enters a two-cycle here:
    [0,0,0] -> [1,1,2] -> [0,1,1] -> [1,1,2]. Both backends counted the cycle
    starting from the step they had just appended — which is the repeat itself,
    so the count stopped at once — and reported every cycle as length 1.
    """
    trace = model.succession([0, 0, 0])
    assert trace["cycle_length"] == 2
    assert trace["converged"] is False


def test_basins_cover_space(model):
    result = model.basins()
    assert result.total == model.n_scenarios
    covered = sum(r["basin_size"] for r in result.scenarios) + result.cycle_count
    assert covered == result.total
    sizes = [r["basin_size"] for r in result.scenarios]
    assert sizes == sorted(sizes, reverse=True)


def test_impact_matrix_shape(model):
    mat = model.impact_matrix()
    n = sum(len(model.variants[d]) for d in model.descriptors)
    assert len(mat) == n and all(len(row) == n for row in mat)


def test_set_impact_edits_in_place(model):
    d1, d2 = model.descriptors[0], model.descriptors[1]
    v1, v2 = model.variants[d1][0], model.variants[d2][1]
    src, tgt = (d1, v1), (d2, v2)

    original = model.get_impact(src, tgt)
    prev = model.set_impact(src, tgt, original + 7)
    assert prev == original
    assert model.get_impact(src, tgt) == original + 7
    assert model.get_impact((0, 0), (1, 1)) == original + 7    # index form, same cell
    model.set_impact(src, tgt, original)
    assert model.get_impact(src, tgt) == original


def test_edit_changes_consistent_set_without_reload(model):
    d1, d2 = model.descriptors[0], model.descriptors[1]
    before = {model.signature(sc) for sc in model.consistent_scenarios()}
    model.set_impact((d2, model.variants[d2][0]), (d1, model.variants[d1][0]), 50)
    model.set_impact((d2, model.variants[d2][-1]), (d1, model.variants[d1][0]), -50)
    after = {model.signature(sc) for sc in model.consistent_scenarios()}
    assert after != before


def test_copy_is_independent(model):
    d1, d2 = model.descriptors[0], model.descriptors[1]
    src, tgt = (d1, model.variants[d1][0]), (d2, model.variants[d2][0])
    baseline = model.get_impact(src, tgt)

    branch = model.copy()
    branch.set_impact(src, tgt, baseline + 100)
    assert branch.get_impact(src, tgt) == baseline + 100
    assert model.get_impact(src, tgt) == baseline     # original untouched


def test_bad_scenario_raises(model):
    with pytest.raises(ValueError):
        model.impact_balance({"nope": "x"})

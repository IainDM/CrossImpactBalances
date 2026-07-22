"""Pure-Python tests for scenario conversion — no Julia runtime required."""

import pytest

from crossimpactbalances import _convert

DESCRIPTORS = ["A", "B", "C"]
VARIANTS = {"A": ["a0", "a1"], "B": ["b0", "b1", "b2"], "C": ["c0", "c1"]}


def test_resolve_name_map():
    u = _convert.resolve_scenario(
        DESCRIPTORS, VARIANTS, {"A": "a1", "B": "b2", "C": "c0"})
    assert u == [1, 2, 0]


def test_resolve_index_map_and_sequence():
    assert _convert.resolve_scenario(
        DESCRIPTORS, VARIANTS, {"A": 0, "B": 1, "C": 1}) == [0, 1, 1]
    assert _convert.resolve_scenario(DESCRIPTORS, VARIANTS, [1, 0, 1]) == [1, 0, 1]


def test_scenario_to_names_roundtrip():
    scenario = {"A": "a1", "B": "b2", "C": "c0"}
    u = _convert.resolve_scenario(DESCRIPTORS, VARIANTS, scenario)
    assert _convert.scenario_to_names(DESCRIPTORS, VARIANTS, u) == scenario


def test_unknown_descriptor():
    with pytest.raises(ValueError, match="Unknown descriptor"):
        _convert.resolve_scenario(DESCRIPTORS, VARIANTS, {"Z": "a0"})


def test_unknown_variant():
    with pytest.raises(ValueError, match="Unknown variant"):
        _convert.resolve_scenario(DESCRIPTORS, VARIANTS, {"A": "nope"})


def test_missing_descriptor():
    with pytest.raises(ValueError, match="Missing descriptor"):
        _convert.resolve_scenario(DESCRIPTORS, VARIANTS, {"A": "a0", "B": "b0"})


def test_out_of_range_index():
    with pytest.raises(ValueError, match="out of range"):
        _convert.resolve_scenario(DESCRIPTORS, VARIANTS, {"A": 5, "B": 0, "C": 0})


def test_wrong_length_sequence():
    with pytest.raises(ValueError, match="descriptors"):
        _convert.resolve_scenario(DESCRIPTORS, VARIANTS, [0, 0])

"""Shared fixtures. Integration tests need the Julia engine; skip gracefully
when it cannot start (e.g. no network to provision Julia)."""

import os

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SAMPLE_DIR = os.path.join(REPO_ROOT, "test", "sample_files")
GLOBAL_SCW = os.path.join(SAMPLE_DIR, "CIB_global.scw")
GLOBAL_SL = os.path.join(SAMPLE_DIR, "CIB_global.sl")


@pytest.fixture(scope="session")
def engine_available():
    """True if the embedded Julia engine can start; otherwise skip."""
    try:
        from crossimpactbalances._engine import get_engine
        get_engine()
        return True
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"Julia engine unavailable: {type(exc).__name__}: {exc}")


@pytest.fixture(scope="session")
def sample_dir():
    return SAMPLE_DIR


@pytest.fixture
def global_model(engine_available):
    from crossimpactbalances import Model
    return Model.load(GLOBAL_SCW, sl_file=GLOBAL_SL)

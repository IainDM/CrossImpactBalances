"""Shared fixtures. Integration tests need an engine backend; skip gracefully
when none can start.

The backend under test is chosen by the ``CIB_BACKEND`` environment variable
(``auto`` by default; ``juliacall`` or ``native`` to pin one) — the same
override :func:`Model.load` consults for ``backend="auto"``, so ``run_models``
and every other entry point agree. The juliacall and native backends each
embed a Julia runtime, so they must not be initialised in the same process —
verify each by running the suite once per backend in separate processes.
"""

import os

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SAMPLE_DIR = os.path.join(REPO_ROOT, "test", "sample_files")
GLOBAL_SCW = os.path.join(SAMPLE_DIR, "CIB_global.scw")

BACKEND = os.environ.get("CIB_BACKEND", "auto")


@pytest.fixture(scope="session")
def engine_available():
    """True if the selected engine backend can load a model; else skip."""
    try:
        from crossimpactbalances import Model
        Model.load(GLOBAL_SCW, backend=BACKEND)
        return True
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"engine backend {BACKEND!r} unavailable: "
                    f"{type(exc).__name__}: {exc}")


@pytest.fixture(scope="session")
def sample_dir():
    return SAMPLE_DIR


@pytest.fixture
def model(engine_available):
    from crossimpactbalances import Model
    return Model.load(GLOBAL_SCW, backend=BACKEND)

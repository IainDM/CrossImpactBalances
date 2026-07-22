"""Python interface to CrossImpactBalances.jl.

Fast Cross-Impact Balance (CIB) scenario analysis, backed by the Julia engine
via `juliacall`. Load a ScenarioWizard ``.scw`` model, find its consistent
scenarios or basins of attraction, score candidate scenarios, and tweak
individual cross-impact values in place without re-parsing.

Quick start
-----------
>>> from crossimpactbalances import Model, run_models, sweep_impact
>>> m = Model.load("model.scw")
>>> m.consistent_scenarios()                      # list of {descriptor: variant}
>>> m.set_impact(("Trade", "Free"), ("Growth", "High"), 3)   # edit in place
>>> m.consistent_scenarios()                      # recomputed, no re-parse
>>> run_models("path/to/scw_dir", analysis="basins")         # batch over files

Multi-threading: set ``PYTHON_JULIACALL_THREADS=auto`` before first use to let
the exhaustive/basin routines run in parallel.
"""

from .batch import run_models, sweep_impact, to_dataframe
from .model import BasinResult, Model

__all__ = [
    "Model",
    "BasinResult",
    "run_models",
    "sweep_impact",
    "to_dataframe",
    "__version__",
]

__version__ = "0.2.0"

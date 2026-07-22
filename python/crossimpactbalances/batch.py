"""Batch helpers — run analyses over many models, or sweep one cross-impact."""

from __future__ import annotations

import glob
import os
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

from .model import Endpoint, Model

_ANALYSES = ("consistent", "basins")


def _expand_paths(paths: Union[str, Iterable[str]]) -> List[str]:
    """Turn a directory, a glob, or an iterable of paths into a list of files."""
    if isinstance(paths, str):
        if os.path.isdir(paths):
            candidates = sorted(glob.glob(os.path.join(paths, "*.scw")))
        else:
            candidates = sorted(glob.glob(paths)) or [paths]
    else:
        candidates = list(paths)
    return candidates


def run_models(paths: Union[str, Iterable[str]], *, analysis: str = "consistent",
               exhaustive: bool = False, algorithm: str = "auto",
               mc_threshold: int = 10000, seed: Optional[int] = None,
               on_error: str = "record", **load_kwargs) -> List[Dict[str, Any]]:
    """Run one analysis over many ``.scw`` model files; one record per file.

    Parameters
    ----------
    paths :
        A directory (its ``*.scw`` files), a glob string, or an iterable of paths.
    analysis :
        ``"consistent"`` (fixed points) or ``"basins"`` (basin analysis).
    on_error :
        ``"record"`` captures a per-file ``{"path", "error"}`` and continues so
        one bad file does not abort the batch; ``"raise"`` re-raises.

    Each record carries the file identity and a structured summary. For
    ``"consistent"``: ``n_consistent`` and the ``scenarios`` list. For
    ``"basins"``: ``cycle_count``, ``largest_basin_fraction`` and the ranked
    ``scenarios``.
    """
    if analysis not in _ANALYSES:
        raise ValueError(f"Unknown analysis {analysis!r}; expected one of {_ANALYSES}")
    if on_error not in ("record", "raise"):
        raise ValueError("on_error must be 'record' or 'raise'")

    records: List[Dict[str, Any]] = []
    for path in _expand_paths(paths):
        try:
            model = Model.load(path, exhaustive=exhaustive,
                               mc_threshold=mc_threshold, seed=seed,
                               compute_kernel=False, **load_kwargs)
            rec: Dict[str, Any] = {
                "path": path,
                "n_descriptors": model.n_descriptors,
                "n_scenarios": model.n_scenarios,
            }
            if analysis == "consistent":
                scenarios = model.consistent_scenarios(
                    exhaustive=exhaustive, algorithm=algorithm, seed=seed)
                rec["n_consistent"] = len(scenarios)
                rec["scenarios"] = scenarios
            else:
                result = model.basins()
                rec["n_consistent"] = result.n_fixed_points
                rec["cycle_count"] = result.cycle_count
                rec["largest_basin_fraction"] = (
                    result.scenarios[0]["basin_fraction"] if result.scenarios else 0.0)
                rec["scenarios"] = result.to_records()
            records.append(rec)
        except Exception as exc:  # noqa: BLE001 - surfaced in the record
            if on_error == "raise":
                raise
            records.append({"path": path, "error": f"{type(exc).__name__}: {exc}"})
    return records


def sweep_impact(model: Model, source: Endpoint, target: Endpoint,
                 values: Sequence[int], *, analysis: str = "consistent",
                 **analysis_kwargs) -> List[Dict[str, Any]]:
    """Vary one cross-impact over *values*, running an analysis for each.

    Edits ``model`` in place for each value (one parse, N scored variants), then
    restores the original value before returning. Each record is
    ``{"value", ...}`` with the analysis summary for that value.
    """
    if analysis not in _ANALYSES:
        raise ValueError(f"Unknown analysis {analysis!r}; expected one of {_ANALYSES}")

    original = model.get_impact(source, target)
    records: List[Dict[str, Any]] = []
    try:
        for value in values:
            model.set_impact(source, target, int(value))
            rec: Dict[str, Any] = {"value": int(value)}
            if analysis == "consistent":
                scenarios = model.consistent_scenarios(**analysis_kwargs)
                rec["n_consistent"] = len(scenarios)
                rec["scenarios"] = scenarios
            else:
                result = model.basins(**analysis_kwargs)
                rec["n_consistent"] = result.n_fixed_points
                rec["cycle_count"] = result.cycle_count
                rec["scenarios"] = result.to_records()
            records.append(rec)
    finally:
        model.set_impact(source, target, original)   # always restore
    return records


def to_dataframe(records: List[Dict[str, Any]]):
    """Flatten batch/sweep records into a pandas ``DataFrame`` (one row per scenario).

    Requires the optional ``pandas`` extra (``pip install
    'crossimpactbalances[pandas]'``). Records that carry an ``"error"`` (or no
    ``"scenarios"``) contribute a single summary row.
    """
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - depends on install extras
        raise ImportError(
            "to_dataframe requires pandas. Install with "
            "'pip install crossimpactbalances[pandas]'.") from exc

    rows: List[Dict[str, Any]] = []
    for rec in records:
        scenarios = rec.get("scenarios")
        base = {k: v for k, v in rec.items() if k != "scenarios"}
        if not scenarios:
            rows.append(base)
            continue
        for i, sc in enumerate(scenarios):
            row = dict(base)
            row["scenario_index"] = i
            # sc is either a {descriptor: variant} map or a basins record.
            row.update(sc)
            rows.append(row)
    return pd.DataFrame(rows)

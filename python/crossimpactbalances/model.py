"""The :class:`Model` — a Pythonic façade over a pluggable engine backend.

Two backends implement the same index-based protocol:

* ``native``    — a compiled ``libcib`` shared library via ctypes; ships no
                  Julia source (see :mod:`._native`).
* ``juliacall`` — the in-process Julia source package (see :mod:`._juliacall`).

:class:`Model` holds whichever backend, caches the model's descriptor/variant
names, and does all human-readable name↔index translation and result
formatting — so callers see the same API regardless of backend.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from . import _convert

Scenario = Union[Mapping[str, Union[str, int]], Sequence[int]]
Endpoint = Tuple[Union[str, int], Union[str, int]]


def _make_backend(path: str, backend: str, load_kwargs: Dict[str, Any]):
    backend = (backend or "auto").lower()
    if backend == "auto":
        # A global override lets callers (and the test suite) pin the backend
        # without threading it through every call — and keeps a single Julia
        # runtime per process.
        backend = os.environ.get("CIB_BACKEND", "auto").lower()
    if backend == "native":
        from ._native import NativeBackend
        return NativeBackend.load(path, **load_kwargs)
    if backend == "juliacall":
        from ._juliacall import JuliaCallBackend
        return JuliaCallBackend.load(path, **load_kwargs)
    if backend == "auto":
        # Prefer the compiled library if one is present; else the source engine.
        try:
            from ._native import _find_library, NativeBackend
            _find_library()
        except Exception:
            from ._juliacall import JuliaCallBackend
            return JuliaCallBackend.load(path, **load_kwargs)
        return NativeBackend.load(path, **load_kwargs)
    raise ValueError(
        f"Unknown backend {backend!r}; expected 'auto', 'native' or 'juliacall'")


class Model:
    """A loaded Cross-Impact Balance model, kept resident by its backend."""

    def __init__(self, backend):
        self._b = backend
        s = backend.structure()
        self._descriptors: List[str] = list(s["descriptors"])
        self._variants: Dict[str, List[str]] = {k: list(v) for k, v in s["variants"].items()}
        self._n_scenarios: int = int(s["n_scenarios"])
        self._path: Optional[str] = None

    # ── construction ────────────────────────────────────────────────────────

    @classmethod
    def load(cls, path: str, *, backend: str = "auto", exhaustive: bool = False,
             mc_threshold: int = 10000, compute_kernel: bool = True,
             sl_file: Optional[str] = None, seed: Optional[int] = None) -> "Model":
        """Parse a ScenarioWizard ``.scw`` file into a :class:`Model`.

        ``backend`` selects the engine: ``"native"`` (compiled library, no
        source), ``"juliacall"`` (in-process Julia source), or ``"auto"`` — the
        native library if one is found, otherwise juliacall.

        ``sl_file`` and ``compute_kernel`` apply to the juliacall backend,
        which can find the kernel at load time; the native backend always
        recomputes on demand and ignores them. ``exhaustive``, ``mc_threshold``
        and ``seed`` are accepted and ignored by both — they configured a
        Monte-Carlo sampling search that no longer exists, the engine now
        always searching exhaustively and deterministically.
        """
        load_kwargs = dict(exhaustive=exhaustive, mc_threshold=mc_threshold,
                           compute_kernel=compute_kernel, sl_file=sl_file, seed=seed)
        m = cls(_make_backend(str(path), backend, load_kwargs))
        m._path = str(path)
        return m

    # ── structure (read-only) ────────────────────────────────────────────────

    @property
    def path(self) -> Optional[str]:
        return self._path

    @property
    def backend(self) -> str:
        return getattr(self._b, "name", "unknown")

    @property
    def descriptors(self) -> List[str]:
        return list(self._descriptors)

    @property
    def variants(self) -> Dict[str, List[str]]:
        return {d: list(v) for d, v in self._variants.items()}

    @property
    def n_descriptors(self) -> int:
        return len(self._descriptors)

    @property
    def n_scenarios(self) -> int:
        return self._n_scenarios

    # ── analysis ──────────────────────────────────────────────────────────────

    def consistent_scenarios(self, *, exhaustive: bool = False,
                             algorithm: str = "auto", ignore_cycles: bool = True,
                             rule: str = "global",
                             seed: Optional[int] = None) -> List[Dict[str, str]]:
        """Consistent scenarios (fixed points) as ``{descriptor: variant}`` maps.

        Recomputes against the current matrix, reflecting any :meth:`set_impact`.
        """
        kern = self._b.find_consistent(exhaustive=exhaustive, algorithm=algorithm,
                                       ignore_cycles=ignore_cycles, rule=rule, seed=seed)
        return [self._names(u) for u in kern]

    def basins(self, *, rule: str = "global") -> "BasinResult":
        """Basin-of-attraction analysis over the whole scenario space."""
        fps, sizes, cycle_count, total = self._b.find_basins(rule=rule)
        rows = []
        for u, size in zip(fps, sizes):
            rows.append({"scenario": self._names(u), "basin_size": size,
                         "basin_fraction": size / total if total else 0.0})
        rows.sort(key=lambda r: r["basin_size"], reverse=True)
        return BasinResult(rows, cycle_count, total)

    def impact_balance(self, scenario: Scenario) -> Dict[str, Any]:
        """Impact scores for *scenario* and whether it is self-consistent."""
        u = self._resolve(scenario)
        ib = self._b.impact_balance(u)
        scores: Dict[str, Dict[str, int]] = {}
        selected: Dict[str, str] = {}
        is_consistent = True
        offset = 0
        for i, desc in enumerate(self._descriptors):
            vars_ = self._variants[desc]
            nv = len(vars_)
            block = ib[offset:offset + nv]
            scores[desc] = {vars_[k]: int(block[k]) for k in range(nv)}
            sel = int(u[i])
            selected[desc] = vars_[sel]
            if block[sel] < max(block):
                is_consistent = False
            offset += nv
        return {"scores": scores, "selected": selected, "is_consistent": is_consistent}

    def is_consistent(self, scenario: Scenario) -> bool:
        """True if *scenario* is a fixed point of the succession dynamics."""
        return self.impact_balance(scenario)["is_consistent"]

    def succession(self, scenario: Scenario, *, rule: str = "global",
                   max_steps: Optional[int] = None) -> Dict[str, Any]:
        """Trace succession from *scenario* to its attractor."""
        u = self._resolve(scenario)
        steps, cycle_length = self._b.succession(u, rule=rule, max_steps=max_steps)
        return {"steps": [self._names(s) for s in steps],
                "converged": cycle_length == 1, "cycle_length": cycle_length}

    def signature(self, scenario: Scenario) -> int:
        """Unique integer signature of *scenario*."""
        return self._b.signature(self._resolve(scenario))

    def scenario_from_signature(self, s: int) -> Dict[str, str]:
        """Inverse of :meth:`signature`: signature -> scenario name map."""
        return self._names(self._b.inv_signature(int(s)))

    # ── in-place editing (no re-parse) ────────────────────────────────────────

    def set_impact(self, source: Endpoint, target: Endpoint, value: int) -> int:
        """Set the cross-impact of *source* variant onto *target* variant.

        ``source``/``target`` are ``(descriptor, variant)`` pairs (name or
        0-based index). Edits the resident model in place and returns the old
        value; the next analysis reflects it — no re-parse.
        """
        sd, sv = self._endpoint(source)
        td, tv = self._endpoint(target)
        return int(self._b.set_impact(sd, sv, td, tv, int(value)))

    def get_impact(self, source: Endpoint, target: Endpoint) -> int:
        """Current cross-impact of *source* variant onto *target* variant."""
        sd, sv = self._endpoint(source)
        td, tv = self._endpoint(target)
        return int(self._b.get_impact(sd, sv, td, tv))

    def impact_matrix(self) -> List[List[int]]:
        """A read-only copy of the full cross-impact matrix. Edit via :meth:`set_impact`."""
        return self._b.matrix()

    def copy(self) -> "Model":
        """An independent deep copy of this model, for trying edits against a baseline."""
        m = Model(self._b.copy())
        m._path = self._path
        return m

    # ── internals ─────────────────────────────────────────────────────────────

    def _names(self, u) -> Dict[str, str]:
        return _convert.scenario_to_names(self._descriptors, self._variants, u)

    def _resolve(self, scenario: Scenario) -> List[int]:
        return _convert.resolve_scenario(self._descriptors, self._variants, scenario)

    def _endpoint(self, pair: Endpoint) -> Tuple[int, int]:
        desc, var = pair
        di = self._desc_index(desc)
        vi = _convert._resolve_variant(self._descriptors[di],
                                       self._variants[self._descriptors[di]], var)
        return di, vi

    def _desc_index(self, desc: Union[str, int]) -> int:
        if isinstance(desc, str):
            try:
                return self._descriptors.index(desc)
            except ValueError:
                raise ValueError(
                    f'Unknown descriptor: "{desc}". '
                    f'Available: {", ".join(self._descriptors)}') from None
        di = int(desc)
        if not (0 <= di < len(self._descriptors)):
            raise ValueError(f"Descriptor index {di} out of range 0:{len(self._descriptors) - 1}")
        return di

    def __repr__(self) -> str:
        return (f"<CIB Model {self._path!r} [{self.backend}]: {self.n_descriptors} "
                f"descriptors, {self.n_scenarios} scenarios>")


class BasinResult:
    """Result of :meth:`Model.basins`: ranked fixed points plus cycle coverage."""

    def __init__(self, scenarios: List[Dict[str, Any]], cycle_count: int, total: int):
        self.scenarios = scenarios          # ranked, largest basin first
        self.cycle_count = cycle_count      # scenarios in cycles (non-convergent)
        self.total = total                  # scenario-space size

    @property
    def n_fixed_points(self) -> int:
        return len(self.scenarios)

    def to_records(self) -> List[Dict[str, Any]]:
        """Flat list of one dict per fixed point (scenario map + basin stats)."""
        out = []
        for rank, row in enumerate(self.scenarios, start=1):
            rec = {"rank": rank, "basin_size": row["basin_size"],
                   "basin_fraction": row["basin_fraction"]}
            rec.update(row["scenario"])
            out.append(rec)
        return out

    def __repr__(self) -> str:
        return (f"<BasinResult: {self.n_fixed_points} fixed points, "
                f"{self.cycle_count} in cycles, {self.total} total>")

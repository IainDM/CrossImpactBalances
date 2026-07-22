"""The :class:`Model` — a Pythonic wrapper around a loaded Julia ``CIB``."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from . import _convert
from ._engine import empty_kernel, get_engine, get_jl, make_rng

# A scenario, Python-side, is a {descriptor: variant} mapping or an index seq.
Scenario = Union[Mapping[str, Union[str, int]], Sequence[int]]
# An endpoint for impact editing: (descriptor, variant), each name or 0-based int.
Endpoint = Tuple[Union[str, int], Union[str, int]]

_RULES = {"global": "GlobalSuccession", "sequential": "SequentialSuccession"}


class Model:
    """A loaded Cross-Impact Balance model.

    Construct via :meth:`Model.load`. The parsed model lives in the embedded
    Julia runtime, so repeated analyses — and in-place edits via
    :meth:`set_impact` — reuse it without re-parsing the ``.scw`` file.
    """

    def __init__(self, cib, path: Optional[str] = None):
        self._cib = cib
        self._path = path
        # Cache the immutable structure as native Python (cheap, converted once).
        self._descriptors: List[str] = [str(d) for d in cib.descriptors]
        self._variants: Dict[str, List[str]] = {
            d: [str(v) for v in cib.variants[d]] for d in self._descriptors
        }

    # ── construction ────────────────────────────────────────────────────────

    @classmethod
    def load(cls, path: str, *, exhaustive: bool = False,
             mc_threshold: int = 10000, compute_kernel: bool = True,
             sl_file: Optional[str] = None, seed: Optional[int] = None) -> "Model":
        """Parse a ScenarioWizard ``.scw`` file into a :class:`Model`.

        Parameters
        ----------
        path : str
            Path to the ``.scw`` model file.
        exhaustive : bool
            Enumerate the full scenario space when computing the initial kernel
            (otherwise a Monte-Carlo sample is used for large spaces).
        mc_threshold : int
            Space size above which Monte-Carlo sampling replaces full enumeration.
        compute_kernel : bool
            If ``False``, skip the load-time kernel search (parse only). The
            kernel is computed on demand by :meth:`consistent_scenarios` anyway.
        sl_file : str, optional
            A ``.sl`` solutions file to seed the kernel from instead of searching.
        seed : int, optional
            Seed for reproducible Monte-Carlo sampling.
        """
        cib_mod = get_engine()
        kwargs: Dict[str, Any] = {
            "mc_threshold": int(mc_threshold),
            "exhaustive": bool(exhaustive),
            "rng": make_rng(seed),
        }
        if sl_file is not None:
            kwargs["sl_file"] = str(sl_file)
        elif not compute_kernel:
            kwargs["kernel"] = empty_kernel()
        cib = cib_mod.load_scw(str(path), **kwargs)
        return cls(cib, path=str(path))

    # ── model structure (read-only) ─────────────────────────────────────────

    @property
    def path(self) -> Optional[str]:
        return self._path

    @property
    def descriptors(self) -> List[str]:
        """Ordered descriptor names."""
        return list(self._descriptors)

    @property
    def variants(self) -> Dict[str, List[str]]:
        """Mapping of descriptor name -> ordered variant names."""
        return {d: list(v) for d, v in self._variants.items()}

    @property
    def n_descriptors(self) -> int:
        return int(self._cib.ndesc)

    @property
    def n_scenarios(self) -> int:
        """Total size of the scenario space (``max_signature + 1``)."""
        return int(get_engine().max_signature(self._cib)) + 1

    # ── analysis ────────────────────────────────────────────────────────────

    def consistent_scenarios(self, *, exhaustive: bool = False,
                             algorithm: str = "auto", ignore_cycles: bool = True,
                             rule: str = "global",
                             seed: Optional[int] = None) -> List[Dict[str, str]]:
        """Return the consistent scenarios (fixed points) as name maps.

        Always recomputes against the current matrix, so it reflects any
        :meth:`set_impact` edits made since loading.
        """
        eng = get_engine()
        kern = eng.find_consistent(
            self._cib,
            rule=self._rule(rule),
            ignore_cycles=bool(ignore_cycles),
            exhaustive=bool(exhaustive),
            algorithm=self._symbol(algorithm),
            rng=make_rng(seed),
        )
        return [self._names(u) for u in kern]

    def basins(self, *, rule: str = "global") -> "BasinResult":
        """Exhaustive basin-of-attraction analysis over the whole scenario space."""
        eng = get_engine()
        fps, sizes, cycle_count = eng.find_basins(self._cib, rule=self._rule(rule))
        total = self.n_scenarios
        rows = []
        for u, size in zip(fps, sizes):
            size = int(size)
            rows.append({
                "scenario": self._names(u),
                "basin_size": size,
                "basin_fraction": size / total if total else 0.0,
            })
        rows.sort(key=lambda r: r["basin_size"], reverse=True)
        return BasinResult(rows, int(cycle_count), total)

    def impact_balance(self, scenario: Scenario) -> Dict[str, Any]:
        """Impact scores for *scenario*, and whether it is self-consistent.

        Returns ``{"scores": {descriptor: {variant: score}}, "selected":
        {descriptor: variant}, "is_consistent": bool}``.
        """
        eng = get_engine()
        u = self._u(scenario)
        ib = _convert.to_index_list(eng.impact_balance(self._cib, u))

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
        """Trace succession from *scenario* to its attractor.

        Returns ``{"steps": [name map, ...], "converged": bool, "cycle_length":
        int}``. ``converged`` is True when the trajectory ends at a fixed point.
        """
        eng = get_engine()
        u = self._u(scenario)
        step_rule = self._rule(rule)
        sig = eng.signature
        start = [int(x) for x in u]
        steps = [start]
        seen = {int(sig(self._cib, u))}
        limit = max_steps if max_steps is not None else self.n_scenarios + 10
        cycle_length = 0
        cur = start
        for _ in range(int(limit)):
            nxt = _convert.to_index_list(
                eng.succession_step(step_rule, self._cib, self._jl_vec(cur)))
            nxt_sig = int(sig(self._cib, self._jl_vec(nxt)))
            if nxt_sig in seen:
                if nxt == steps[-1]:
                    cycle_length = 1
                else:
                    steps.append(nxt)
                    for k in range(len(steps) - 1, -1, -1):
                        cycle_length += 1
                        if int(sig(self._cib, self._jl_vec(steps[k]))) == nxt_sig:
                            break
                break
            seen.add(nxt_sig)
            steps.append(nxt)
            cur = nxt
        return {
            "steps": [self._names(s) for s in steps],
            "converged": cycle_length == 1,
            "cycle_length": cycle_length,
        }

    def signature(self, scenario: Scenario) -> int:
        """Unique integer signature of *scenario*."""
        return int(get_engine().signature(self._cib, self._u(scenario)))

    def scenario_from_signature(self, s: int) -> Dict[str, str]:
        """Inverse of :meth:`signature`: signature -> scenario name map."""
        u = get_engine().inv_signature(self._cib, int(s))
        return self._names(u)

    # ── in-place editing (no re-parse) ──────────────────────────────────────

    def set_impact(self, source: Endpoint, target: Endpoint, value: int) -> int:
        """Set the cross-impact of *source* variant onto *target* variant.

        ``source`` and ``target`` are ``(descriptor, variant)`` pairs, each part
        given by name or 0-based index. Edits the resident model in place (both
        the matrix and its stored transpose) and returns the previous value.
        The next :meth:`consistent_scenarios` / :meth:`basins` call reflects the
        change — no re-parse of the ``.scw``.
        """
        sd, sv = source
        td, tv = target
        old = get_engine().set_impact_b(
            self._cib, self._norm(sd), self._norm(sv),
            self._norm(td), self._norm(tv), int(value))
        return int(old)

    def get_impact(self, source: Endpoint, target: Endpoint) -> int:
        """Current cross-impact of *source* variant onto *target* variant."""
        sd, sv = source
        td, tv = target
        return int(get_engine().get_impact(
            self._cib, self._norm(sd), self._norm(sv),
            self._norm(td), self._norm(tv)))

    def impact_matrix(self) -> List[List[int]]:
        """A read-only copy of the full cross-impact matrix (``cim``).

        Rows/columns are flat variant indices in descriptor order. Edit via
        :meth:`set_impact`, not by mutating this copy.
        """
        cim = self._cib.cim
        return [[int(x) for x in row] for row in cim]

    def copy(self) -> "Model":
        """A deep copy of this model, for trying edits against a baseline."""
        return Model(get_jl().deepcopy(self._cib), path=self._path)

    # ── internals ───────────────────────────────────────────────────────────

    def _names(self, u) -> Dict[str, str]:
        return _convert.scenario_to_names(self._descriptors, self._variants, u)

    def _u(self, scenario: Scenario):
        """Resolve a Python scenario to a Julia ``Vector{Int}``."""
        idx = _convert.resolve_scenario(self._descriptors, self._variants, scenario)
        return _convert.to_julia_index_vector(idx)

    def _jl_vec(self, indices: Sequence[int]):
        return _convert.to_julia_index_vector(indices)

    def _rule(self, rule: str):
        try:
            ctor = _RULES[rule]
        except KeyError:
            raise ValueError(
                f"Unknown rule {rule!r}; expected one of {sorted(_RULES)}") from None
        return get_jl().seval(f"CrossImpactBalances.{ctor}()")

    @staticmethod
    def _symbol(name: str):
        if name not in ("auto", "bnb", "sweep"):
            raise ValueError(
                f"Unknown algorithm {name!r}; expected 'auto', 'bnb' or 'sweep'")
        return get_jl().seval(f":{name}")

    @staticmethod
    def _norm(part: Union[str, int]):
        """Normalise a descriptor/variant token for the Julia helpers."""
        return part if isinstance(part, str) else int(part)

    def __repr__(self) -> str:
        return (f"<CIB Model {self._path!r}: {self.n_descriptors} descriptors, "
                f"{self.n_scenarios} scenarios>")


class BasinResult:
    """Result of :meth:`Model.basins`: ranked fixed points plus cycle coverage."""

    def __init__(self, scenarios: List[Dict[str, Any]], cycle_count: int,
                 total: int):
        self.scenarios = scenarios          # ranked, largest basin first
        self.cycle_count = cycle_count      # scenarios in cycles (non-convergent)
        self.total = total                  # scenario-space size

    @property
    def n_fixed_points(self) -> int:
        return len(self.scenarios)

    def to_records(self) -> List[Dict[str, Any]]:
        """Flat list of one dict per fixed point (scenario dict + basin stats)."""
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

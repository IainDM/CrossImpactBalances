"""juliacall backend — drive the engine in-process via the Julia source package.

Implements the same index-based protocol as :class:`NativeBackend`, so
:class:`Model` uses either interchangeably. This backend needs the Julia
source package available (it develops ``CrossImpactBalances.jl``); use it for
development in this repository. To ship without exposing source, use the
native (compiled-library) backend instead.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import _convert
from ._engine import empty_kernel, get_engine, get_jl, make_rng

_RULES = {"global": "GlobalSuccession", "sequential": "SequentialSuccession"}


class JuliaCallBackend:
    name = "juliacall"

    def __init__(self, cib):
        self._cib = cib
        eng = get_engine()
        descriptors = [str(d) for d in cib.descriptors]
        self._structure = {
            "descriptors": descriptors,
            "variants": {d: [str(v) for v in cib.variants[d]] for d in descriptors},
            "n_descriptors": int(cib.ndesc),
            "n_scenarios": int(eng.max_signature(cib)) + 1,
        }

    @classmethod
    def load(cls, path: str, *, exhaustive: bool = False, mc_threshold: int = 10000,
             compute_kernel: bool = True, sl_file: Optional[str] = None,
             seed: Optional[int] = None) -> "JuliaCallBackend":
        eng = get_engine()
        kwargs: Dict[str, Any] = {
            "mc_threshold": int(mc_threshold),
            "exhaustive": bool(exhaustive),
            "rng": make_rng(seed),
        }
        if sl_file is not None:
            kwargs["sl_file"] = str(sl_file)
        elif not compute_kernel:
            kwargs["kernel"] = empty_kernel()
        return cls(eng.load_scw(str(path), **kwargs))

    def structure(self) -> Dict[str, Any]:
        return self._structure

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

    def _vec(self, indices):
        return _convert.to_julia_index_vector(indices)

    def find_consistent(self, *, exhaustive=False, algorithm="auto",
                        ignore_cycles=True, rule="global", seed=None) -> List[List[int]]:
        kern = get_engine().find_consistent(
            self._cib, rule=self._rule(rule), ignore_cycles=bool(ignore_cycles),
            exhaustive=bool(exhaustive), algorithm=self._symbol(algorithm),
            rng=make_rng(seed))
        return [[int(x) for x in u] for u in kern]

    def find_basins(self, *, rule="global"):
        eng = get_engine()
        fps, sizes, cycle_count = eng.find_basins(self._cib, rule=self._rule(rule))
        total = int(eng.max_signature(self._cib)) + 1
        return ([[int(x) for x in u] for u in fps],
                [int(s) for s in sizes], int(cycle_count), total)

    def impact_balance(self, u: List[int]) -> List[int]:
        return [int(x) for x in get_engine().impact_balance(self._cib, self._vec(u))]

    def succession(self, u: List[int], *, rule="global", max_steps=None):
        eng = get_engine()
        sig = eng.signature
        step_rule = self._rule(rule)
        start = [int(x) for x in u]
        steps = [start]
        seen = {int(sig(self._cib, self._vec(start)))}
        limit = max_steps if max_steps is not None else self._structure["n_scenarios"] + 10
        cycle_length = 0
        cur = start
        for _ in range(int(limit)):
            nxt = [int(x) for x in eng.succession_step(step_rule, self._cib, self._vec(cur))]
            nsig = int(sig(self._cib, self._vec(nxt)))
            if nsig in seen:
                if nxt == steps[-1]:
                    cycle_length = 1
                else:
                    steps.append(nxt)
                    for k in range(len(steps) - 1, -1, -1):
                        cycle_length += 1
                        if int(sig(self._cib, self._vec(steps[k]))) == nsig:
                            break
                break
            seen.add(nsig)
            steps.append(nxt)
            cur = nxt
        return steps, cycle_length

    def signature(self, u: List[int]) -> int:
        return int(get_engine().signature(self._cib, self._vec(u)))

    def inv_signature(self, s: int) -> List[int]:
        return [int(x) for x in get_engine().inv_signature(self._cib, int(s))]

    def set_impact(self, sd, sv, td, tv, value) -> int:
        return int(get_engine().set_impact_b(
            self._cib, int(sd), int(sv), int(td), int(tv), int(value)))

    def get_impact(self, sd, sv, td, tv) -> int:
        return int(get_engine().get_impact(
            self._cib, int(sd), int(sv), int(td), int(tv)))

    def matrix(self) -> List[List[int]]:
        # Iterating a Julia Matrix yields scalars column-major; build rows in
        # Julia (as the C-API does) so both backends return the same shape.
        rows = get_jl().seval("c -> [c.cim[i, :] for i in 1:size(c.cim, 1)]")(self._cib)
        return [[int(x) for x in row] for row in rows]

    def copy(self) -> "JuliaCallBackend":
        return JuliaCallBackend(get_jl().deepcopy(self._cib))

    def close(self):
        pass

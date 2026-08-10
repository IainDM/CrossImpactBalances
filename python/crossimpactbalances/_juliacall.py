"""juliacall backend — drive the engine in-process via the Julia source package.

Implements the same index-based protocol as :class:`NativeBackend`, so
:class:`Model` uses either interchangeably. This backend needs the Julia
source package available (it develops ``CrossImpactBalances.jl``); use it for
development in this repository. To ship without exposing source, use the
native (compiled-library) backend instead.

Counts and signatures are read out of Julia through their decimal digits (see
:func:`_exact_int`), because a real ScenarioWizard model can have more
scenarios than ``Int64`` can count — 10^24 is an ordinary size — and both
backends have to agree about numbers that large. The C API crosses the same
values as decimal strings for the same reason; Python's ``int`` is
arbitrary-precision, so here they simply arrive intact.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import _convert
from ._engine import empty_kernel, get_engine, get_jl

_RULES = {"global": "GlobalSuccession", "sequential": "SequentialSuccession"}

_INT64_MAX = 2 ** 63 - 1


def _exact_int(value) -> int:
    """A Julia integer as a Python ``int``, exactly, whatever its width.

    ``scenario_count`` and ``_signature128`` return ``Int128``. Going via the
    decimal digits keeps every one of them, where anything routed through a
    float would round past 2**53 — and real models reach 10**30.
    """
    return int(str(value))


class JuliaCallBackend:
    name = "juliacall"

    def __init__(self, cib):
        self._cib = cib
        eng = get_engine()
        descriptors = [str(d) for d in cib.descriptors]
        self._structure = {
            "descriptors": descriptors,
            "variants": {d: [str(v) for v in cib.variants[d]] for d in descriptors},
            "n_descriptors": int(cib.numberOfDescriptors),
            # scenario_count, not max_signature + 1: max_signature THROWS past
            # typemax(Int64), which would make a model that size unloadable
            # here even though find_consistent searches it perfectly well.
            "n_scenarios": _exact_int(eng.scenario_count(cib)),
        }

    @classmethod
    def load(cls, path: str, *, compute_kernel: bool = True,
             sl_file: Optional[str] = None, **_ignored) -> "JuliaCallBackend":
        """Parse *path*, mirroring :meth:`NativeBackend.load`.

        ``_ignored`` absorbs the load options :class:`Model` still accepts but
        the engine no longer has — ``exhaustive``, ``mc_threshold`` and
        ``seed`` configured a Monte-Carlo sampling search that was removed when
        the search became exhaustive and deterministic. The native backend
        drops them the same way.
        """
        eng = get_engine()
        kwargs: Dict[str, Any] = {}
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

    def _sig128(self, u) -> int:
        """The Int128 signature of *u*, as an exact Python int.

        ``CrossImpactBalances.signature`` returns an ``Int`` and wraps silently
        into negative numbers past ``typemax(Int64)``; ``_signature128`` is the
        same mixed-radix number in a type wide enough for any loadable model.
        It is Julia-private, hence the underscore — the engine keys its own
        large-space tallies on it.
        """
        return _exact_int(get_engine()._signature128(self._cib, self._vec(u)))

    def find_consistent(self, *, algorithm="auto", rule="global",
                        bnb_node_budget=None, **_ignored) -> List[List[int]]:
        """The model's kernel. ``_ignored`` absorbs the removed search options
        (``exhaustive``, ``ignore_cycles``, ``seed``), as the native backend does."""
        kwargs: Dict[str, Any] = {"rule": self._rule(rule),
                                  "algorithm": self._symbol(algorithm)}
        if bnb_node_budget is not None:
            kwargs["bnb_node_budget"] = int(bnb_node_budget)
        kern = get_engine().find_consistent(self._cib, **kwargs)
        return [[int(x) for x in u] for u in kern]

    def find_basins(self, *, rule="global"):
        eng = get_engine()
        fps, sizes, cycle_count = eng.find_basins(self._cib, rule=self._rule(rule))
        total = _exact_int(eng.scenario_count(self._cib))
        return ([[int(x) for x in u] for u in fps],
                [int(s) for s in sizes], int(cycle_count), total)

    def impact_balance(self, u: List[int]) -> List[int]:
        return [int(x) for x in get_engine().impact_balance(self._cib, self._vec(u))]

    def succession(self, u: List[int], *, rule="global", max_steps=None):
        eng = get_engine()
        step_rule = self._rule(rule)
        start = [int(x) for x in u]
        steps = [start]
        # Int128 signatures as the "have I been here before" key: past
        # typemax(Int64) the Int64 ones wrap, two different scenarios can share
        # a key, and the walk reports a cycle that is not there.
        seen = {self._sig128(start)}
        limit = max_steps if max_steps is not None else self._structure["n_scenarios"] + 10
        cycle_length = 0
        cur = start
        for _ in range(int(limit)):
            nxt = [int(x) for x in eng.succession_step(step_rule, self._cib, self._vec(cur))]
            nsig = self._sig128(nxt)
            if nsig in seen:
                if nxt == steps[-1]:
                    cycle_length = 1
                else:
                    steps.append(nxt)
                    # Count back from the step BEFORE the one just appended:
                    # starting AT it matches immediately (nsig is its own
                    # signature), so every cycle came back as length 1 — which
                    # Model.succession turns into ``converged: True``.
                    for k in range(len(steps) - 2, -1, -1):
                        cycle_length += 1
                        if self._sig128(steps[k]) == nsig:
                            break
                break
            seen.add(nsig)
            steps.append(nxt)
            cur = nxt
        return steps, cycle_length

    def signature(self, u: List[int]) -> int:
        return self._sig128(u)

    def inv_signature(self, s: int) -> List[int]:
        s = int(s)
        if s > _INT64_MAX:
            # Same refusal the C API gives, for the same reason: signatures are
            # reported exactly at any size, but inverting one is Int64 work.
            raise ValueError(
                f"signature {s} is past typemax(Int64) = {_INT64_MAX}, and "
                "inv_signature works in Int64 — it cannot invert it. The "
                "scenario's variant indices identify it at any model size.")
        return [int(x) for x in get_engine().inv_signature(self._cib, s)]

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

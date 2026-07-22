"""Boundary conversions between Julia values and native Python types.

Scenarios are represented on the Python side as ``{descriptor: variant}`` name
maps. Internally the Julia engine uses length-``ndesc`` vectors of 0-based
variant indices. These helpers translate between the two, mirroring the
0-based/1-based conventions in ``src/CrossImpactBalances.jl`` and the
name-resolution logic in ``mcp/julia_worker.jl``.
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Sequence

from ._engine import get_jl


def to_index_list(jl_vector) -> List[int]:
    """Convert a Julia ``Vector{Int}`` to a Python ``list[int]``."""
    return [int(x) for x in jl_vector]


def scenario_to_names(descriptors: Sequence[str],
                      variants: Mapping[str, Sequence[str]],
                      u) -> Dict[str, str]:
    """Index vector ``u`` (0-based) -> ``{descriptor: variant_name}`` dict.

    Mirrors ``cib.variants[desc][u[i] + 1]`` (Julia 1-based) as ``[u[i]]`` in
    Python's 0-based indexing.
    """
    out: Dict[str, str] = {}
    for i, desc in enumerate(descriptors):
        out[desc] = variants[desc][int(u[i])]
    return out


def resolve_scenario(descriptors: Sequence[str],
                     variants: Mapping[str, Sequence[str]],
                     scenario) -> List[int]:
    """Resolve a scenario to a length-``ndesc`` list of 0-based variant indices.

    *scenario* may be:
      - a ``{descriptor: variant}`` mapping (variant given by name, or by int
        index), OR
      - a sequence of 0-based indices already in descriptor order.

    Every descriptor must be specified. Raises ``ValueError`` on unknown
    descriptor/variant names, out-of-range indices, or missing descriptors —
    with the same "Available: ..." hints as the Julia/worker code.
    """
    ndesc = len(descriptors)

    # Sequence form: already in descriptor order.
    if not isinstance(scenario, Mapping):
        seq = list(scenario)
        if len(seq) != ndesc:
            raise ValueError(
                f"Scenario has {len(seq)} entries but the model has {ndesc} "
                f"descriptors: {list(descriptors)}")
        out = []
        for i, desc in enumerate(descriptors):
            out.append(_resolve_variant(desc, variants[desc], seq[i]))
        return out

    # Mapping form: {descriptor: variant}.
    index = [-1] * ndesc
    desc_pos = {d: i for i, d in enumerate(descriptors)}
    for key, val in scenario.items():
        if key not in desc_pos:
            raise ValueError(
                f'Unknown descriptor: "{key}". '
                f'Available: {", ".join(descriptors)}')
        i = desc_pos[key]
        index[i] = _resolve_variant(key, variants[key], val)

    missing = [descriptors[i] for i in range(ndesc) if index[i] == -1]
    if missing:
        raise ValueError(
            "Missing descriptor(s): " + ", ".join(f'"{m}"' for m in missing) +
            ". All descriptors must be specified.")
    return index


def _resolve_variant(desc: str, vars_: Sequence[str], value) -> int:
    """Resolve a variant given by name or 0-based index to a 0-based index."""
    if isinstance(value, str):
        try:
            return list(vars_).index(value)
        except ValueError:
            raise ValueError(
                f'Unknown variant "{value}" for descriptor "{desc}". '
                f'Available: {", ".join(vars_)}') from None
    v = int(value)
    if not (0 <= v < len(vars_)):
        raise ValueError(
            f"Variant index {v} out of range 0:{len(vars_) - 1} for "
            f'descriptor "{desc}"')
    return v


def to_julia_index_vector(indices: Sequence[int]):
    """Build a Julia ``Vector{Int}`` from a Python sequence of ints."""
    jl = get_jl()
    return jl.seval("Vector{Int}")([int(i) for i in indices])

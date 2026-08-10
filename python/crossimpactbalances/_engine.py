"""Lazy Julia runtime bootstrap for the CrossImpactBalances Python interface.

Importing :mod:`juliacall` starts an embedded Julia, and (via juliapkg, using
the ``juliapkg.json`` shipped next to this module) provisions Julia itself and
develops the local ``CrossImpactBalances.jl`` package. That is expensive, so it
is deferred until the first real use rather than paid at ``import
crossimpactbalances`` time.

Multi-threading: the fast exhaustive/basin paths only parallelise when Julia is
started with several threads. Set the environment variable
``PYTHON_JULIACALL_THREADS=auto`` (or a specific count) *before* the first call
into the engine to enable it.
"""

from __future__ import annotations

import os
import threading

_LOCK = threading.Lock()
_JL = None      # the juliacall Main module
_CIB = None     # the CrossImpactBalances Julia module handle


def _prepare_env():
    """Set juliacall env vars that must be in place before it is imported.

    When Julia is started with multiple threads, juliacall requires
    ``PYTHON_JULIACALL_HANDLE_SIGNALS=yes`` — otherwise multithreaded routines
    (find_basins, exhaustive search) can segfault, especially when another
    signal handler (e.g. pytest, a debugger) is installed. We opt in
    automatically whenever more than one thread is requested, unless the user
    has already set the variable. (The trade-off: Julia then handles signals,
    so Python's Ctrl-C won't raise KeyboardInterrupt.)
    """
    threads = os.environ.get("PYTHON_JULIACALL_THREADS", "").strip().lower()
    wants_threads = threads not in ("", "1")
    if wants_threads and "PYTHON_JULIACALL_HANDLE_SIGNALS" not in os.environ:
        os.environ["PYTHON_JULIACALL_HANDLE_SIGNALS"] = "yes"


def _bootstrap():
    global _JL, _CIB
    _prepare_env()
    from juliacall import Main as jl  # heavy: starts Julia, resolves deps

    jl.seval("using CrossImpactBalances")
    _JL = jl
    _CIB = jl.CrossImpactBalances


def get_jl():
    """Return the juliacall ``Main`` module, starting Julia on first call."""
    if _JL is None:
        with _LOCK:
            if _JL is None:
                _bootstrap()
    return _JL


def get_engine():
    """Return the ``CrossImpactBalances`` Julia module, starting Julia if needed."""
    if _CIB is None:
        with _LOCK:
            if _CIB is None:
                _bootstrap()
    return _CIB


def empty_kernel():
    """A concretely-typed empty ``Vector{Vector{Int}}`` for ``load_scw(kernel=...)``.

    A bare Python ``[]`` converts to a Julia ``Vector{Any}`` and would fail
    dispatch, so build the concrete type explicitly.
    """
    return get_jl().seval("Vector{Vector{Int}}()")

"""
    CrossImpactBalances

Cross-Impact Balance (CIB) scenario analysis. A Julia implementation of the methodology in Weimer-Jehle (2006), ported from the Python [sei-international/cibsa](https://github.com/sei-international/cibsa) library with a two-phase basin-of-attraction analysis, a threaded incremental exhaustive sweep, and an exact branch-and-bound search that prunes provably-inconsistent regions of the scenario space. The succession dynamics can be changed through [`SuccessionRule`](@ref). There is a default called [`GlobalSuccession`](@ref) which follows the standard ScenarioWizard rule.

# Reference
Weimer-Jehle, W. (2006). Cross-impact balances: A system-theoretical approach to cross-impact analysis. *Technological Forecasting and Social Change*, 73(4), 334–361.
"""
module CrossImpactBalances

# `export` lists the names visible after `using CrossImpactBalances` — like 
# Python's `__all__` or C#'s `public`. Everything else in the package is still
# reachable as `CrossImpactBalances.name`, it just isn't brought into scope
# automatically. Names starting with an underscore are internal helpers.
export CIB, load_scw, load_solutions,
       impact_balance,
       set_impact!, get_impact,
       SuccessionRule, GlobalSuccession, SequentialSuccession,
       succession_step,
       find_consistent, find_basins,
       signature, inv_signature, max_signature

# ══ Vocabulary ═════════════════════════════════════════════════════════════
#
# Three terms in this codebase mean the same thing. They are kept because each
# comes from a different place and readers arrive expecting different ones:
#
#   CONSISTENT SCENARIO  The CIB term (Weimer-Jehle). A scenario in which no
#                        descriptor would rather switch to a different variant.
#                        This is the DEFAULT term used in prose throughout.
#
#   FIXED POINT          The dynamical-systems term for the same thing: a
#                        scenario whose succession step maps it to itself.
#                        Used where that framing carries weight — chiefly in
#                        the basin analysis, where the meaningful contrast is
#                        fixed point versus cycle.
#
#   KERNEL               The collective noun: the set of all consistent
#                        scenarios in a model. It appears in the public API
#                        (`load_scw(...; kernel=...)`) so it is kept.
#
# So "find the kernel", "find every consistent scenario" and "find every fixed
# point of the succession map" are three ways of saying one thing.
#
# ── Julia notes for readers coming from Python or .NET ─────────────────────
#
# Why Julia? Because it's fast, in a nutshell. In tests we got a ~30x speedup vs Python
# before any algorithmic speedups.

# Features that recur throughout these files:
#
# * 1-based indexing: Julia arrays start at index 1 (like MATLAB/Fortran),
#   and a range `1:n` includes BOTH endpoints. Scenarios themselves store
#   0-based variant numbers (matching the Python original and the signature
#   arithmetic), which is why many array accesses add `+ 1`.
#
# * Multiple dispatch: several functions share one name but differ in their
#   argument TYPES (e.g. `succession_step(::GlobalSuccession, ...)` vs
#   `succession_step(::SequentialSuccession, ...)`). Julia picks the method
#   from the runtime types of all arguments — like C# method overloading,
#   but resolved dynamically, and it is the standard way to write what
#   Python/C# would express with an interface or virtual method.
#
# * A trailing `!` in a function name (e.g. `push!`, `set_impact!`) is a
#   naming convention meaning "this function modifies one of its arguments".
#   It has no effect on behaviour — it is a warning label.
#
# * Words starting with `@` are macros — code transformations applied before
#   compilation. The ones used here are performance/threading annotations:
#   `@inbounds` (skip array bounds checks), `@simd` (allow vectorisation),
#   `Threads.@spawn` (run on a worker thread, like Task.Run), `@sync` (wait
#   for all tasks spawned inside the block, like Task.WaitAll), and `@view`
#   (slice an array without copying it).
#
# * `:auto`, `:sweep`, `:bnb` are Symbols — lightweight interned names, used
#   here the way C# would use an enum or Python a short string constant.
#
# * `Union{Nothing, Int}` is a nullable type (C#'s `int?`); `nothing` plays
#   the role of null/None, tested with `isnothing(x)` or `x === nothing`.
#
# * `condition && action` and `condition || action` are used as one-line if
#   statements: `&&` runs the action only when the condition is true,
#   `||` only when it is false (e.g. `isempty(s) && continue`).
#
# * `"""docstrings"""` directly above a definition are attached to it as
#   documentation (visible via `?name` in the REPL); `[`Name`](@ref)` inside
#   them is a cross-reference link for the documentation generator.

# ══ The files, in reading order ════════════════════════════════════════════
#
# `include` splices a file in as though its text were typed here, so all of
# these share this one module and its exports — they are NOT separate Julia
# modules and need no imports of each other. Order matters only in that a
# type must be defined before it is named in another file's signatures.

# ── The model, and how to get one ──
include("types.jl")            # the CIB struct: what a model IS
include("signatures.jl")       # numbering every scenario 0, 1, 2, ...
include("impact_balance.jl")   # scoring every variant against a scenario
include("parsers.jl")          # reading ScenarioWizard .scw / .sl files
include("matrix_editing.jl")   # reading/editing one impact, no re-parse

# ── The dynamics: how one scenario leads to the next ──
include("succession.jl")       # SuccessionRule and the built-in rules

# ── Analysis 1: which scenarios are consistent? ──
include("find_consistent.jl")  # the entry point, and choosing a strategy
include("sweep.jl")            #   strategy A: check every scenario, fast
include("branch_and_bound.jl") #   strategy B: prune whole families at once

# ── Analysis 2: where does each scenario end up? ──
include("basins.jl")           # basins of attraction

end # module

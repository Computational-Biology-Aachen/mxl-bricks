# ADR 0001: One File Per Reaction "Brick", Composed Through a Canonical `names` Registry

**Status:** Implemented
**Scope:** `src/mxlbricks/enzymes/*.py`, `src/mxlbricks/names.py`, `src/mxlbricks/fns.py`

---

## 1. Context

`mxlbricks` builds mechanistic ODE models by composing dozens of independent reaction
"bricks" (`enzymes/atp_synthase.py`, `enzymes/rubisco_poolman.py`, ...) into full models
in `models.py`. Every brick needs to refer to shared species/parameter names (e.g. `ATP`,
`e0` for an enzyme's total concentration) consistently — a typo'd string literal for a
species name in one brick silently fails to connect to the same species used by another
brick, since `mxlpy` models key everything by string name.

## 2. Decision

**2.1 — One reaction per file under `enzymes/`.** Each brick is an `add_*(model, ...) ->
Model` function in its own module, taking a `Model` and mutating it (adding
variables/parameters/reactions), following the same fluent, `Self`-returning convention
as `mxlpy`'s own `Model` (see
[mxlpy ADR 0002](https://github.com/Computational-Biology-Aachen/MxlPy/blob/main/docs/adrs/0002-fluent-builder-and-cache-invalidation.md)).

**2.2 — All names come from `names.py`, never inline string literals.** Species,
parameter, and enzyme names are constructed via typed functions in `names.py` (e.g.
`n.atp()`, `n.e0(n.rubisco())`), not hand-typed strings scattered across brick files.
Shared rate-law functions live in `fns.py` and are named, not written as inline lambdas
(inheriting [mxlpy ADR 0004](https://github.com/Computational-Biology-Aachen/MxlPy/blob/main/docs/adrs/0004-named-functions-no-lambdas.md)).

## 3. Rationale

The one-brick-per-file layout makes each reaction independently reviewable, testable, and
reusable across models — `models.py`'s `get_*` functions are essentially a recipe of
which bricks to combine, which is only legible if each brick is a small, self-contained
unit. Routing all names through `names.py` turns "two bricks refer to the same species"
from a stringly-typed hope into something a type checker and IDE autocomplete can verify
— a typo in `n.atp()` is a `NameError`/import error, not a silently-disconnected species
that only surfaces as a wrong simulation result.

## 4. Consequences

- Adding a new brick means adding both the `add_*` function under `enzymes/` *and* any
  new name constructors it needs to `names.py` — skipping the latter and inlining a
  string literal defeats the whole point and should be treated as a bug, not a shortcut.
- Because every brick takes and returns the same `model: Model` object, bricks compose
  freely in any order (subject to their own dependencies), which is what makes
  `models.py`'s dense `get_*` functions ([ADR 0004](0004-drifted-relationship-with-mxlmodels.md))
  tractable to read despite wiring together 20+ bricks.

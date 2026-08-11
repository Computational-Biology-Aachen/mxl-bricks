# ADR 0003: `enzymes/__init__.py` Is Generated, Not Hand-Maintained

**Status:** Implemented
**Scope:** `src/mxlbricks/enzymes/_create_imports.py`, `src/mxlbricks/enzymes/__init__.py`

---

## 1. Context

`enzymes/` holds 80+ brick modules, each exporting one or more public `add_*` functions
that `models.py` needs to import. Hand-maintaining an `__init__.py` re-exporting all of
them as bricks are added, renamed, or removed is exactly the kind of bookkeeping that
silently drifts out of sync (a new brick file added without updating `__init__.py`
compiles fine but is simply unreachable from `models.py`).

## 2. Decision

`_create_imports.py` is a small standalone script (stdlib `ast`, no dependency on
`mxlbricks` itself) that walks every non-private module in `enzymes/`, extracts its
top-level public function names, and rewrites `__init__.py` from scratch as
`from .<module> import (<fn1>, <fn2>, ...)` lines plus a computed `__all__`. It is run
manually (not on every build) after adding or changing bricks.

## 3. Rationale

Generating `__init__.py` mechanically from the actual module contents makes "did I
forget to export the new brick" structurally impossible to get wrong — the script is the
single source of truth for what's exported, derived from the code itself rather than
maintained in parallel with it. Keeping it a small `ast`-based script with no framework
dependency (rather than, say, a `mxlpy` metaprogramming feature) matches the narrow scope
of the problem: this is a code-organization convenience specific to `mxlbricks`'
`enzymes/` layout, not part of `mxlpy`'s own `meta/` codegen system
(see [mxlpy ADR 0008](https://github.com/Computational-Biology-Aachen/MxlPy/blob/main/docs/adrs/0008-meta-codegen-single-source-of-truth.md)),
which projects a *model* to multiple targets — this instead projects a *file layout* to
an import list.

## 4. Consequences

- Never hand-edit `enzymes/__init__.py` directly — the next run of
  `_create_imports.py` will silently overwrite manual changes. Add/rename/remove a brick
  module, then re-run the script.
- Because the script is not wired into CI or a pre-commit hook, a forgotten re-run after
  adding a brick will fail loudly at import time in `models.py` (the new `add_*` function
  won't be importable) rather than silently — acceptable given how infrequently new
  bricks are added, but worth automating if that cadence increases.

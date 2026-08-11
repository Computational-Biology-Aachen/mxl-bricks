# mxlbricks: Architecture Context

This is the entry point for understanding *why* `mxlbricks` is shaped the way it is —
written down ahead of a maintainer handoff, alongside the equivalent `docs/adrs/`
directories in the sibling `mxlpy`, `mxlmodels`, `absorpig`, `parameteriser`, and
`schemegen` repos.

## Composing Models from Bricks

`mxlbricks` builds mechanistic ODE models by assembling standalone, reusable reaction
"bricks" via `mxlpy`'s `Model` builder — the whole point is sharing reaction machinery
across a lineage of related photosynthesis models rather than re-deriving it per model.

→ [ADR 0001 — One file per reaction brick, composed through a canonical `names` registry](0001-brick-per-file-with-canonical-name-registry.md)
→ [ADR 0003 — `enzymes/__init__.py` is generated, not hand-maintained](0003-codegen-enzymes-init.md)

## Trustworthy Performance Optimizations

→ [ADR 0002 — `matrix` / `analytical` / `analytical-split`: three cross-validated formulations of the same kinetics](0002-matrix-analytical-analytical-split-variants.md)

## Relationship to `mxlmodels`

→ [ADR 0004 — Relationship to `mxlmodels` has drifted: codegen source for some models, not all](0004-drifted-relationship-with-mxlmodels.md)

This is the fact most likely to surprise a new contributor: `mxlbricks` is *not* the
source of everything in `mxlmodels`. Only the shared photosynthesis-core lineage
(Yokota 1985 → Poolman 2000 → Matuszyńska 2016/2019 → Saadat 2021 → Ebeling 2026) is
actually assembled from bricks here; `mxlmodels` has grown well beyond that lineage on
its own.

## Threads That Cross Multiple ADRs

- **Mechanical correctness over hand-derived cleverness, verified not assumed.**
  ADR 0002's matrix-as-ground-truth stance and ADR 0003's generated `__init__.py` are the
  same idea applied twice: prefer a form that's mechanically derived (from equations, or
  from the file system) over a hand-maintained one, and when a hand-crafted shortcut
  exists anyway (a closed-form kinetic expression), keep a permanent automated check that
  it still agrees with the mechanical source of truth.
- **Decomposition into bricks is a cost, not a default.** ADR 0004's drifted relationship
  with `mxlmodels` reflects that brick decomposition only pays off for models that
  actually share machinery with siblings — don't expect or impose it on every model.

## See Also

- [`mxlpy`'s `docs/adrs/CONTEXT.md`](https://github.com/Computational-Biology-Aachen/MxlPy/blob/main/docs/adrs/CONTEXT.md),
  especially ADR 0002 (fluent builder) and ADR 0004 (named functions, no lambdas) —
  conventions `mxlbricks`' bricks directly inherit.
- [`mxlmodels`' `docs/adrs/`](https://github.com/Computational-Biology-Aachen/mxl-models/tree/main/docs/adrs)
  for the same brick/flat-file relationship from the consuming side, and for the many
  models that live only there.

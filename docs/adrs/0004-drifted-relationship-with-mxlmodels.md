# ADR 0004: Relationship to `mxlmodels` Has Drifted — Codegen Source for Some Models, Not All

**Status:** Accepted (documenting an evolved reality, not a new decision)
**Scope:** cross-repo — `mxlbricks/src/mxlbricks/__init__.py` vs.
`mxlmodels/src/mxlmodels/__init__.py`

---

## 1. Context

`mxlbricks` exports exactly seven full-model builders:
`get_ebeling_2026`, `get_matuszynska2016_npq`, `get_matuszynska2016_phd`,
`get_matuszynska2019`, `get_poolman2000`, `get_saadat2021`, `get_yokota1985`. `mxlmodels`
ships single-file, flat versions of models for easy inspection, and historically these
flat files were produced by codegen *from* `mxlbricks`' brick composition. Today,
`mxlmodels` contains over 25 models, most of which (`bellasio2019`, `davis2017`,
`lotka_volterra_*`, `sir`, `prigogine1968_brusselator`, `elowitz2000_repressilator`,
`hahn1987`, `lazar1997`, `zhu2005`, `zhu2009`, the `ss/` steady-state models, ...) have no
corresponding brick composition in `mxlbricks` at all — they were written directly as
flat `mxlpy` models.

## 2. Decision

Accept that `mxlbricks` → `mxlmodels` codegen is the origin story for only the
photosynthesis-core lineage that both packages export under matching names (Yokota 1985,
Poolman 2000, the Matuszyńska 2016/2019 family, Saadat 2021, and now Ebeling 2026) — not
a general rule that every `mxlmodels` model has a `mxlbricks` counterpart. Do not assume
structural symmetry between the two packages' contents.

## 3. Rationale

`mxlbricks`' brick-composition approach earns its complexity when a model shares
substantial reaction machinery with siblings in the same lineage (each new photosynthesis
model in this family reuses most of the previous one's bricks, adding only a few new
ones) — that's precisely the set of models it actually covers. Models added to
`mxlmodels` outside that lineage (classic ODE benchmarks like Lotka-Volterra/SIR/
Brusselator, or one-off published models with no shared brick vocabulary) have nothing to
gain from decomposition into reusable bricks, so they were written directly as flat
`mxlpy` models in `mxlmodels` instead of being retrofitted into `mxlbricks`.

## 4. Consequences

- When asked to "find the `mxlbricks` version" of an arbitrary `mxlmodels` model, check
  `mxlbricks/src/mxlbricks/__init__.py`'s export list first — most `mxlmodels` models
  simply don't have one, and that's expected, not a gap to fill.
- Adding a new model to the shared photosynthesis lineage should still go through
  `mxlbricks` bricks first, with the `mxlmodels` flat file generated from it — but a
  standalone/one-off model added to `mxlmodels` has no obligation to be decomposed into
  `mxlbricks` bricks.
- See `mxlmodels`' own ADR on this same relationship
  ([`docs/adrs/0001-flat-files-mixed-provenance.md`](https://github.com/Computational-Biology-Aachen/mxl-models/blob/main/docs/adrs/0001-flat-files-mixed-provenance.md))
  for the inverse view — which models are codegen'd vs. hand-written, from the
  `mxlmodels` side.

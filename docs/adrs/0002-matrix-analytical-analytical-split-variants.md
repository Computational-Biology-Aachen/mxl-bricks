# ADR 0002: `matrix` / `analytical` / `analytical-split` — Three Cross-Validated Formulations of the Same Kinetics

**Status:** Implemented
**Scope:** `_m16_npq.py`, `models.py` (`get_matuszynska2016_npq`, `get_matuszynska2016_phd`,
`get_matuszynska2019`, `get_saadat2021`), `tests/test_variants.py`

---

## 1. Context

Several photosystem-state models (2016 NPQ, 2016 PhD, 2019, Saadat 2021) involve fast
quasi-steady-state photosystem-state kinetics. These can be expressed either as a linear
system solved numerically each step (`mxlpy.surrogates.qss`, the `"matrix"` mode — the
form derived directly from the conservation/mass-action equations) or as a hand-derived
closed-form algebraic expression (`"analytical"`/`"analytical-split"` modes) that is
faster to evaluate and easier to read, at the cost of being far more error-prone to
derive by hand.

## 2. Decision

Every such model exposes a `mode: Literal["matrix", "analytical", "analytical-split"]`
parameter that selects between the three formulations, and
`tests/test_variants.py` asserts numerically (`atol=1e-12, rtol=1e-12`) that all three
produce identical `get_args()` output for the same model — including at high light
(`pfd=5000`), where the closed forms are most likely to diverge from the matrix ground
truth if the hand-derivation has an error.

## 3. Rationale

The `"matrix"` form is treated as ground truth because it is mechanically derived from
the QSSA linear system rather than hand-manipulated algebra, so it is the form least
likely to contain a derivation error — but it's also the slowest to evaluate at
simulation time. The closed-form variants exist purely for performance and readability,
and are only trustworthy once verified to agree with the matrix form; keeping the test
suite's numeric cross-check permanently in place (not just as a one-time derivation
check) guards against a future refactor of either form silently breaking the equivalence.
The dedicated high-light test exists because these are exactly the extreme-input
conditions where a subtly-wrong closed form is most likely to diverge from the
mechanically-correct matrix form.

## 4. Consequences

- Never hand-edit an `"analytical"`/`"analytical-split"` variant without re-running
  `test_variants.py` — a passing test is the only thing establishing that the closed form
  is actually equivalent to the matrix ground truth, not just plausible-looking algebra.
- When adding a new model with this kind of fast sub-system, default to implementing the
  `"matrix"` form first and treat closed-form variants as an optional, separately-verified
  optimization — not the other way around.
- `"analytical"` vs `"analytical-split"` is itself a further trade-off (e.g. one grouped
  expression vs. split intermediate terms) — consult the specific model file for which
  is preferred for readability vs. reuse in that context.

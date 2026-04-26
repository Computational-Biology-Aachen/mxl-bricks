"""phosphohexomutase.

EC 5.3.1.9

Equilibrator
D-Fructose 6-phosphate(aq) ⇌ D-Glucose 6-phosphate(aq)
Keq = 3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_1s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_glucose_6_phosphate_isomerase_re(
    model: Model,
    *,
    rxn: str | None = None,
    f6p: str | None = None,
    g6p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glucose 6 phosphate isomerase re (phosphohexomutase) to model."""
    rxn = default_name(rxn, n.g6pi)
    f6p = default_name(f6p, n.f6p)
    g6p = default_name(g6p, n.g6p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_1s_1p,
        stoichiometry={
            f6p: -1,
            g6p: 1,
        },
        args=[
            f6p,
            g6p,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=2.3),
        ],
    )
    return model

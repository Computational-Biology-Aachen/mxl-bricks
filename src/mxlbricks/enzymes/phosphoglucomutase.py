"""glucose phosphomutase.

EC 5.4.2.2

G6P <=> G1P

Equilibrator
Glucose 6-phosphate(aq) ⇌ D-Glucose-1-phosphate(aq)
Keq = 0.05 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_1s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_phosphoglucomutase(
    model: Model,
    *,
    rxn: str | None = None,
    g6p: str | None = None,
    g1p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add phosphoglucomutase (glucose phosphomutase) to model."""
    rxn = default_name(rxn, n.phosphoglucomutase)
    g6p = default_name(g6p, n.g6p)
    g1p = default_name(g1p, n.g1p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_1s_1p,
        stoichiometry={
            g6p: -1,
            g1p: 1,
        },
        args=[
            g6p,
            g1p,
            default_kre(model, rxn=rxn, par=kre, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=0.058),
        ],
    )
    return model

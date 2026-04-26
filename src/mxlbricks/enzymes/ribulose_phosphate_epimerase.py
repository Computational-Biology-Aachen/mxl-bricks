"""ribulose-phosphate 3-epimerase.

EC 5.1.3.1

Equilibrator
    D-Xylulose 5-phosphate(aq) ⇌ D-Ribulose 5-phosphate(aq)
    Keq = 0.3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_1s_1p
from mxlbricks.utils import (
    default_kre,
    default_name,
    static,
)


def add_ribulose_5_phosphate_3_epimerase(
    model: Model,
    *,
    rxn: str | None = None,
    x5p: str | None = None,
    ru5p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add ribulose 5 phosphate 3 epimerase (ribulose-phosphate 3-epimerase) to model."""
    rxn = default_name(rxn, n.ribulose_phosphate_epimerase)
    x5p = default_name(x5p, n.x5p)
    ru5p = default_name(ru5p, n.ru5p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_1s_1p,
        stoichiometry={
            x5p: -1,
            ru5p: 1,
        },
        args=[
            x5p,
            ru5p,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            static(model, n.keq(rxn), 0.67) if keq is None else keq,
        ],
    )
    return model

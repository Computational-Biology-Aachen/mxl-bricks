"""ribose-5-phosphate isomerase.

EC 5.3.1.6

Equilibrator
    D-Ribose 5-phosphate(aq) ⇌ D-Ribulose 5-phosphate(aq)
    Keq = 0.4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_1s_1p
from mxlbricks.utils import (
    default_kre,
    default_name,
    static,
)


def add_ribose_5_phosphate_isomerase(
    model: Model,
    *,
    rxn: str | None = None,
    r5p: str | None = None,
    ru5p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add ribose 5 phosphate isomerase (ribose-5-phosphate isomerase) to model."""
    rxn = default_name(rxn, n.ribose_phosphate_isomerase)
    r5p = default_name(r5p, n.r5p)
    ru5p = default_name(ru5p, n.ru5p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_1s_1p,
        stoichiometry={
            r5p: -1,
            ru5p: 1,
        },
        args=[
            r5p,
            ru5p,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            static(model, n.keq(rxn), 0.4) if keq is None else keq,
        ],
    )
    return model

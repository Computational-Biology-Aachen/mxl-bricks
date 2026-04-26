"""triose-phosphate isomerase.

EC 5.3.1.1

Equilibrator
    D-Glyceraldehyde 3-phosphate(aq) ⇌ Glycerone phosphate(aq)
    Keq = 10 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_1s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_triose_phosphate_isomerase(
    model: Model,
    *,
    rxn: str | None = None,
    gap: str | None = None,
    dhap: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add triose phosphate isomerase (triose-phosphate isomerase) to model."""
    rxn = default_name(rxn, n.triose_phosphate_isomerase)
    gap = default_name(gap, n.gap)
    dhap = default_name(dhap, n.dhap)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_1s_1p,
        stoichiometry={
            gap: -1,
            dhap: 1,
        },
        args=[
            gap,
            dhap,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=22.0),
        ],
    )
    return model

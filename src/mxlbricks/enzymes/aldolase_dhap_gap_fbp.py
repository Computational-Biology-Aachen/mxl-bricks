"""DHAP + GAP <=> FBP.

EC 4.1.2.13

Equilibrator
Glycerone phosphate(aq) + D-Glyceraldehyde 3-phosphate(aq) ⇌ D-Fructose 1,6-bisphosphate(aq)
Keq = 1.1e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)

"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_2s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_aldolase_dhap_gap_req(
    model: Model,
    *,
    rxn: str | None = None,
    gap: str | None = None,
    dhap: str | None = None,
    fbp: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add aldolase dhap gap req (DHAP + GAP <=> FBP) to model."""
    rxn = default_name(rxn, n.aldolase_dhap_gap)
    gap = default_name(gap, n.gap)
    dhap = default_name(dhap, n.dhap)
    fbp = default_name(fbp, n.fbp)
    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_1p,
        stoichiometry={
            gap: -1,
            dhap: -1,
            fbp: 1,
        },
        args=[
            gap,
            dhap,
            fbp,
            default_kre(model, rxn=rxn, par=kre, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=7.1),
        ],
    )
    return model

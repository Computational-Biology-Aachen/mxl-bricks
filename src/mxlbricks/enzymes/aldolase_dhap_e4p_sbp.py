"""DHAP + EAP <=> SBP.

EC 4.1.2.13

Equilibrator
Glycerone phosphate(aq) + D-Erythrose 4-phosphate(aq) ⇌ Sedoheptulose 1,7-bisphosphate(aq)
Keq = 4.8e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_2s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_aldolase_dhap_e4p_req(
    model: Model,
    *,
    rxn: str | None = None,
    dhap: str | None = None,
    e4p: str | None = None,
    sbp: str | None = None,
    keq: str | None = None,
    kre: str | None = None,
) -> Model:
    """Add aldolase dhap e4p req (DHAP + EAP <=> SBP) to model."""
    rxn = default_name(rxn, n.aldolase_dhap_e4p)
    dhap = default_name(dhap, n.dhap)
    e4p = default_name(e4p, n.e4p)
    sbp = default_name(sbp, n.sbp)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_1p,
        stoichiometry={
            dhap: -1,
            e4p: -1,
            sbp: 1,
        },
        args=[
            dhap,
            e4p,
            sbp,
            default_kre(model, rxn=rxn, par=kre, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=13.0),
        ],
    )
    return model

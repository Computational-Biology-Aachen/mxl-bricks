"""DHAP + A5P <=> OBP.

Equilibrator
metacyc.compound:DIHYDROXY-ACETONE-PHOSPHATE + metacyc.compound:ARABINOSE-5P
    ⇌ metacyc.compound:CPD-17017
Keq = 0.0509 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_2s_1p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_aldolase_dhap_a5p_req(
    model: Model,
    *,
    rxn: str | None = None,
    dhap: str | None = None,
    a5p: str | None = None,
    obp: str | None = None,
    keq: str | None = None,
    kre: str | None = None,
) -> Model:
    """Add aldolase dhap a5p req (DHAP + A5P <=> OBP) to model."""
    rxn = default_name(rxn, n.aldolase_dhap_e4p)
    dhap = default_name(dhap, n.dhap)
    a5p = default_name(a5p, n.arabinose_5_phosphate)
    obp = default_name(obp, n.obp)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_1p,
        stoichiometry={
            dhap: -1,
            a5p: -1,
            obp: 1,
        },
        args=[
            dhap,
            a5p,
            obp,
            default_kre(model, rxn=rxn, par=kre, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=0.0509),
        ],
    )
    return model

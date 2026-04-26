"""name.

EC 6.4.1.2

Metacyc:
ACETYL-COA_m + ATP_m + HCO3_m <=> ADP_m + MALONYL-COA_m + PROTON_m + Pi_m

Equilibrator
Acetyl-CoA(aq) + ATP(aq) + HCO3-(aq) ⇌ ADP(aq) + Malonyl-CoA(aq) + Orthophosphate(aq)
Too much uncertainty for HCO3

As a proxy
Acetyl-CoA(aq) + ATP(aq) + CO2(total) ⇌ ADP(aq) + Malonyl-CoA(aq) + Orthophosphate(aq)
Keq = 4e1 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    reversible_michaelis_menten_3s_3p,
    reversible_michaelis_menten_3s_3p_1i,
)
from mxlbricks.utils import (
    default_keq,
    default_ki,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_acetyl_coa_carboxyltransfer(
    model: Model,
    *,
    rxn: str | None = None,
    s1: str | None = None,
    s2: str | None = None,
    s3: str | None = None,
    p1: str | None = None,
    p2: str | None = None,
    p3: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add acetyl coa carboxyltransfer (name) to model."""
    rxn = default_name(rxn, n.acetyl_coa_carboxyltransfer)
    s1 = default_name(s1, n.acetyl_coa)
    s2 = default_name(s2, n.atp)
    s3 = default_name(s3, n.hco3)
    p1 = default_name(p1, n.adp)
    p2 = default_name(p2, n.malonyl_coa)
    p3 = default_name(p3, n.pi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            stoichiometry={
                s1: -1.0,
                s2: -1.0,
                s3: -1.0,
                p1: 1.0,
                p2: 1.0,
                p3: 1.0,
            },
        ),
        args=[
            s1,
            s2,
            s3,
            p1,
            p2,
            p3,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=30.1,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.0487),
            default_kmp(model, rxn=rxn, par=kmp, value=0.1),
            default_keq(model, rxn=rxn, par=keq, value=40.0),
        ],
    )

    return model


def add_acetyl_coa_carboxyltransfer_1i(
    model: Model,
    *,
    rxn: str | None = None,
    s1: str | None = None,
    s2: str | None = None,
    s3: str | None = None,
    p1: str | None = None,
    p2: str | None = None,
    p3: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
    i1: str | None = None,
    ki: str | None = None,
) -> Model:
    """Add acetyl coa carboxyltransfer 1i (name) to model."""
    rxn = default_name(rxn, n.acetyl_coa_carboxyltransfer)
    s1 = default_name(s1, n.acetyl_coa)
    s2 = default_name(s2, n.atp)
    s3 = default_name(s3, n.hco3)
    p1 = default_name(p1, n.adp)
    p2 = default_name(p2, n.malonyl_coa)
    p3 = default_name(p3, n.pi)
    i1 = default_name(i1, n.formate)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p_1i,
        stoichiometry=filter_stoichiometry(
            model,
            stoichiometry={
                s1: -1.0,
                s2: -1.0,
                s3: -1.0,
                p1: 1.0,
                p2: 1.0,
                p3: 1.0,
            },
        ),
        args=[
            s1,
            s2,
            s3,
            p1,
            p2,
            p3,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=30.1,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.0487),
            default_kmp(model, rxn=rxn, par=kmp, value=0.1),
            default_keq(model, rxn=rxn, par=keq, value=40.0),
            default_ki(model, rxn=rxn, par=ki, value=0.002),
        ],
    )

    return model

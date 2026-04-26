"""name.

EC 2.3.1.9

Metacyc:
ACETOACETYL-COA_m + CO-A_m  <=>  2.0 ACETYL-COA_m

Equilibrator
Acetoacetyl-CoA(aq) + CoA(aq) ⇌ 2 Acetyl-CoA(aq)
Keq = 2.4e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_2s_1p
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_acetyl_coa_acetyltransfer(
    model: Model,
    *,
    rxn: str | None = None,
    acac: str | None = None,
    coa: str | None = None,
    ac_coa: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add acetyl coa acetyltransfer (name) to model."""
    rxn = default_name(rxn, n.acetyl_coa_acetyltransfer)
    acac = default_name(acac, n.acetoacetyl_coa)
    coa = default_name(coa, n.coa)
    ac_coa = default_name(ac_coa, n.acetyl_coa)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            stoichiometry={
                acac: -1,
                coa: -1,
                ac_coa: 2,
            },
        ),
        args=[
            acac,
            coa,
            ac_coa,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=220.5,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.0176),
            default_kmp(model, rxn=rxn, par=kmp, value=0.1386),
            default_keq(model, rxn=rxn, par=keq, value=24000.0),
        ],
    )

    return model

"""EC 1.2.1.75.

Metacyc:
MALONYL-COA_m + NADPH_m + PROTON_m <=> CO-A_m + MALONATE-S-ALD_m + NADP_m

Equilibrator
Malonyl-CoA(aq) + NADPH(aq) ⇌ Malonate semialdehyde(aq) + NADP(aq) + CoA(aq)
Keq = 5.6e-3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_2s_3p
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


def add_malonyl_coa_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    malonyl_coa: str | None = None,
    nadph: str | None = None,
    malonate_s_aldehyde: str | None = None,
    nadp: str | None = None,
    coa: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add malonyl coa reductase (EC 1) to model."""
    rxn = default_name(rxn, n.malonyl_coa_reductase)
    malonyl_coa = default_name(malonyl_coa, n.malonyl_coa)
    nadph = default_name(nadph, n.nadph)
    malonate_s_aldehyde = default_name(malonate_s_aldehyde, n.malonate_s_aldehyde)
    nadp = default_name(nadp, n.nadp)
    coa = default_name(coa, n.coa)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                malonyl_coa: -1,
                nadph: -1,
                malonate_s_aldehyde: 1,
                coa: 1,
            },
        ),
        args=[
            malonyl_coa,
            nadph,
            malonate_s_aldehyde,
            nadp,
            coa,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=50.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.03),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=0.0056),
        ],
    )
    return model

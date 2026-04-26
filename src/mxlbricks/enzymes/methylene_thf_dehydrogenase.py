"""EC 1.5.1.5.

Metacyc: METHYLENETHFDEHYDROG-NADP-RXN
METHENYL-THF_m + NADPH_m + 0.93 PROTON_m <=> METHYLENE-THF_m + NADP_m

Equilibrator
5,10-Methenyltetrahydrofolate(aq) + NADPH(aq) ⇌ 5,10-Methylenetetrahydrofolate(aq) + NADP(aq)
Keq = 1e1 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_2s_2p
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


def add_methylene_thf_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    methenyl_thf: str | None = None,
    nadph: str | None = None,
    methylene_thf: str | None = None,
    nadp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add methylene thf dehydrogenase (EC 1) to model."""
    rxn = default_name(rxn, n.methylene_thf_dehydrogenase)
    methenyl_thf = default_name(methenyl_thf, n.methenyl_thf)
    nadph = default_name(nadph, n.nadph)
    methylene_thf = default_name(methylene_thf, n.methylene_thf)
    nadp = default_name(nadp, n.nadp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                methenyl_thf: -1,
                nadph: -1,
                methylene_thf: 1,
                nadp: 1,
            },
        ),
        args=[
            methenyl_thf,
            nadph,
            methylene_thf,
            nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=14.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.12),
            default_kmp(model, rxn=rxn, par=kmp, value=0.302),
            default_keq(model, rxn=rxn, par=keq, value=10.0),
        ],
    )

    return model

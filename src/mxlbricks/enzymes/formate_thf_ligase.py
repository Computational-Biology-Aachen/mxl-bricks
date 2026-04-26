"""name.

EC 6.3.4.3
    Metacyc: FORMATETHFLIG-RXN
    FORMATE_m + ATP_m + THF_m <=> 10-FORMYL-THF_m + ADP_m + Pi_m

    Equilibrator
    Formate(aq) + THF(aq) + ATP(aq) ⇌ 10-Formyltetrahydrofolate(aq) + ADP(aq) + Pi(aq)
    Keq = 2.0 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_3s_3p
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


def add_formate_thf_ligase(
    model: Model,
    *,
    rxn: str | None = None,
    formate: str | None = None,
    atp: str | None = None,
    thf: str | None = None,
    formyl_thf: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add formate thf ligase (name) to model."""
    rxn = default_name(rxn, n.formate_thf_ligase)
    formate = default_name(formate, n.formate)
    atp = default_name(atp, n.atp)
    thf = default_name(thf, n.thf)
    formyl_thf = default_name(formyl_thf, n.formyl_thf)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                formate: -1,
                atp: -1,
                thf: -1,
                formyl_thf: 1,
                adp: 1,
                pi: 1,
            },
        ),
        args=[
            formate,
            atp,
            thf,
            formyl_thf,
            adp,
            pi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=6.08,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=7.6),
            default_kmp(model, rxn=rxn, par=kmp, value=10.0),
            default_keq(model, rxn=rxn, par=keq, value=2.0),
        ],
    )

    return model

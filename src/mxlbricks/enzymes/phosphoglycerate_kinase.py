"""Phosphoglycerate kinase (PGK).

EC 2.7.2.3

kcat
    - 537 | 1 /s | Pseudomonas sp. | brenda

km
    - 0.18 | PGA | mM | Synechocystis sp. | brenda
    - ? | BPGA | mM | Synechocystis sp. | brenda
    - 0.3 | ATP | mM | Spinacia oleracea | brenda
    - 0.27 | ADP | mM | Spinacia oleracea | brenda


Equilibrator
    ATP(aq) + 3-Phospho-D-glycerate(aq) ⇌ ADP(aq) + 3-Phospho-D-glyceroyl phosphate(aq)
    Keq = 3.7e-4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    rapid_equilibrium_2s_2p,
    reversible_michaelis_menten_2s_2p,
)
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_kre,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_phosphoglycerate_kinase_poolman(
    model: Model,
    *,
    rxn: str | None = None,
    pga: str | None = None,
    atp: str | None = None,
    bpga: str | None = None,
    adp: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add phosphoglycerate kinase poolman (Phosphoglycerate kinase (PGK)) to model."""
    rxn = default_name(rxn, n.phosphoglycerate_kinase)
    pga = default_name(pga, n.pga)
    atp = default_name(atp, n.atp)
    bpga = default_name(bpga, n.bpga)
    adp = default_name(adp, n.adp)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                pga: -1.0,
                atp: -1.0,
                bpga: 1.0,
                adp: 1.0,
            },
        ),
        args=[
            pga,
            atp,
            bpga,
            adp,
            default_kre(model, rxn=rxn, par=kre, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=0.00031),
        ],
    )
    return model


def add_phosphoglycerate_kinase(
    model: Model,
    *,
    rxn: str | None = None,
    pga: str | None = None,
    atp: str | None = None,
    bpga: str | None = None,
    adp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add phosphoglycerate kinase (Phosphoglycerate kinase (PGK)) to model."""
    rxn = default_name(rxn, n.phosphoglycerate_kinase)
    pga = default_name(pga, n.pga)
    atp = default_name(atp, n.atp)
    bpga = default_name(bpga, n.bpga)
    adp = default_name(adp, n.adp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                pga: -1.0,
                atp: -1.0,
                bpga: 1.0,
                adp: 1.0,
            },
        ),
        args=[
            pga,
            atp,
            bpga,
            adp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=537,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.18),
            default_kmp(model, rxn=rxn, par=kmp, value=0.27),
            default_keq(model, rxn=rxn, par=keq, value=3.7e-4),
        ],
    )
    return model

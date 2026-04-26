"""glycerate kinase.

ATP + D-Glycerate <=> ADP + 3-Phospho-D-glycerate

Equilibrator
ATP(aq) + D-Glycerate(aq) ⇌ ADP(aq) + 3-Phospho-D-glycerate(aq)
Keq = 4.9e2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_glycerate_kinase(
    model: Model,
    *,
    rxn: str | None = None,
    glycerate: str | None = None,
    atp: str | None = None,
    pga: str | None = None,
    adp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycerate kinase (glycerate kinase) to model."""
    rxn = default_name(rxn, n.glycerate_kinase)
    glycerate = default_name(glycerate, n.glycerate)
    atp = default_name(atp, n.atp)
    pga = default_name(pga, n.pga)
    adp = default_name(adp, n.adp)

    # FIXME: km_atp missing
    # FIXME: ki missing
    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry={
            glycerate: -1.0,
            atp: -1.0,
            pga: 1.0,
        },
        args=[
            glycerate,
            atp,
            pga,
            adp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=5.71579,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.25),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=490.0),
        ],
    )
    return model

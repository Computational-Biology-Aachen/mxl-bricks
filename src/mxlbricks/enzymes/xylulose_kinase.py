"""XYLULOSE + ATP  <=> XYLULOSE_5_PHOSPHATE + ADP.

EC FIXME

Equilibrator
Xylulose(aq) + ATP(aq) ⇌ Xylulose-5-phosphate(aq) + ADP(aq)
Keq = 2.2e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_xylulose_kinase(
    model: Model,
    *,
    rxn: str | None = None,
    xylulose: str | None = None,
    atp: str | None = None,
    x5p: str | None = None,
    adp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add xylulose kinase (XYLULOSE + ATP  <=> XYLULOSE_5_PHOSPHATE + ADP) to model."""
    rxn = default_name(rxn, n.xylulose_kinase)
    xylulose = default_name(xylulose, n.xylulose)
    atp = default_name(atp, n.atp)
    x5p = default_name(x5p, n.x5p)
    adp = default_name(adp, n.adp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                xylulose: -1.0,
                atp: -1.0,
                x5p: 1.0,
                adp: 1.0,
            },
        ),
        args=[
            xylulose,
            atp,
            x5p,
            adp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.1),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=2.2e4),
        ],
    )
    return model

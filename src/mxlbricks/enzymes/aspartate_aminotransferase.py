"""Aspartate aminotransferase.

EC 2.6.1.1

Equilibrator
Aspartate(aq) + alpha-Ketoglutarate(aq) ⇌ Oxaloacetate(aq) + Glutamate(aq)
Keq = 0.3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_aspartate_aminotransferase(
    model: Model,
    *,
    rxn: str | None = None,
    s1: str | None = None,
    s2: str | None = None,
    p1: str | None = None,
    p2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add aspartate aminotransferase (Aspartate aminotransferase) to model."""
    rxn = default_name(rxn, n.aspartate_aminotransferase)
    s1 = default_name(s1, n.aspartate)
    s2 = default_name(s2, n.oxoglutarate)
    p1 = default_name(p1, n.oxaloacetate)
    p2 = default_name(p2, n.glutamate)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                s1: -1.0,
                s2: -1.0,
                p1: 1.0,
                p2: 1.0,
            },
        ),
        args=[
            s1,
            s2,
            p1,
            p2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=84,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=2.53),
            default_kmp(model, rxn=rxn, par=kmp, value=3.88),
            default_keq(model, rxn=rxn, par=keq, value=0.3),
        ],
    )
    return model

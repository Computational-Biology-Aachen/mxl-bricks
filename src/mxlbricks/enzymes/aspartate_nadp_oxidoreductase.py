"""Aspartate:NADP oxidoreductase.

EC FIXME

Equilibrator
Iminoaspartate(aq) + NADPH(aq) ⇌ Aspartate(aq) + NADP(aq)
Keq = 1.6e10 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_2s
from mxlbricks.utils import (
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_aspartate_nadp_oxidoreductase(
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
) -> Model:
    """Add aspartate nadp oxidoreductase (Aspartate:NADP oxidoreductase) to model."""
    rxn = default_name(rxn, n.aspartate_oxidoreductase)
    s1 = default_name(s1, n.iminoaspartate)
    s2 = default_name(s2, n.nadph)
    p1 = default_name(p1, n.aspartate)
    p2 = default_name(p2, n.nadp)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
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
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.1),
        ],
    )
    return model

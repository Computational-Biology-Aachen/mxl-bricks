"""name.

EC FIXME

Equilibrator
Glycolaldehyde(aq) + DHAP(aq) ⇌ Ribose-1-phosphate(aq)
Keq = 0.5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_r1p_aldolase(
    model: Model,
    rxn: str | None = None,
    s1: str | None = None,
    s2: str | None = None,
    p1: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add r1p aldolase (name) to model."""
    rxn = default_name(rxn, n.r1p_aldolase)
    s1 = default_name(s1, n.glycolaldehyde)
    s2 = default_name(s2, n.dhap)
    p1 = default_name(p1, n.r1p)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                s1: -1.0,
                s2: -1.0,
                p1: 1.0,
            },
        ),
        args=[
            s1,
            s2,
            p1,
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
            default_keq(model, rxn=rxn, par=keq, value=0.5),
        ],
    )
    return model

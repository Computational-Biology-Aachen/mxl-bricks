"""Spontaneous reaction.

EC FIXME

Equilibrator
Iminoaspartate(aq) + H2O(l) ⇌ Oxaloacetate(aq) + NH3(aq)
Keq = 7.1e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_1s_2p
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


def add_oxaloacetate_formation(
    model: Model,
    *,
    rxn: str | None = None,
    iminoaspartate: str | None = None,
    oxaloacetate: str | None = None,
    nh4: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add oxaloacetate formation (Spontaneous reaction) to model."""
    rxn = default_name(rxn, n.oxaloacetate_formation)
    iminoaspartate = default_name(iminoaspartate, n.iminoaspartate)
    oxaloacetate = default_name(oxaloacetate, n.oxaloacetate)
    nh4 = default_name(nh4, n.nh4)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                iminoaspartate: -1,
                oxaloacetate: 1,
                nh4: 1,
            },
        ),
        args=[
            iminoaspartate,
            oxaloacetate,
            nh4,
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
            default_keq(model, rxn=rxn, par=keq, value=7100.0),
        ],
    )
    return model

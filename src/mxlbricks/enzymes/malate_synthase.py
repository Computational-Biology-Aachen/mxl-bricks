"""malate synthase.

EC 2.3.3.9

Equilibrator
------------
    H2O + Acetyl-CoA + Glyoxylate <=> CoA + (S)-Malate
    Keq = 6.0e6 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)

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


def add_malate_synthase(
    model: Model,
    *,
    rxn: str | None = None,
    acetyl_coa: str | None = None,
    glyoxylate: str | None = None,
    coa: str | None = None,
    malate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add malate synthase (malate synthase) to model."""
    rxn = default_name(rxn, n.malate_synthase)
    acetyl_coa = default_name(acetyl_coa, n.acetyl_coa)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    coa = default_name(coa, n.coa)
    malate = default_name(malate, n.malate)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                acetyl_coa: -1,
                glyoxylate: -1,
                coa: 1,
                malate: 1,
            },
        ),
        args=[
            acetyl_coa,
            glyoxylate,
            coa,
            malate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=27.8,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.098),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=6.0e6),
        ],
    )
    return model

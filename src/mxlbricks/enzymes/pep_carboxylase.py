"""name.

EC 4.1.1.31

Equilibrator
Phosphoenolpyruvate(aq) + CO2(total) ⇌ Orthophosphate(aq) + Oxaloacetate(aq)
Keq = 4.4e5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_pep_carboxylase(
    model: Model,
    *,
    rxn: str | None = None,
    pep: str | None = None,
    hco3: str | None = None,
    oxaloacetate: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add pep carboxylase (name) to model."""
    rxn = default_name(rxn, n.pep_carboxylase)
    pep = default_name(pep, n.pep)
    hco3 = default_name(hco3, n.hco3)
    oxaloacetate = default_name(oxaloacetate, n.oxaloacetate)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            stoichiometry={
                pep: -1,
                hco3: -1,
                pi: 1,
                oxaloacetate: 1,
            },
        ),
        args=[
            pep,
            hco3,
            oxaloacetate,
            pi,
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
            default_keq(model, rxn=rxn, par=keq, value=440000.0),
        ],
    )
    return model

"""name.

EC 1.1.1.37

Equilibrator
Oxaloacetate(aq) + NADH(aq) ⇌ Malate(aq) + NAD(aq)
Keq = 4.4e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_malate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    oxaloacetate: str | None = None,
    nadh: str | None = None,
    malate: str | None = None,
    nad: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add malate dehydrogenase (name) to model."""
    rxn = default_name(rxn, n.malate_dehydrogenase)
    oxaloacetate = default_name(oxaloacetate, n.oxaloacetate)
    nadh = default_name(nadh, n.nadh)
    malate = default_name(malate, n.malate)
    nad = default_name(nad, n.nad)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                oxaloacetate: -1,
                nadh: -1,
                malate: 1,
                nad: 1,
            },
        ),
        args=[
            oxaloacetate,
            nadh,
            malate,
            nad,
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
            default_keq(model, rxn=rxn, par=keq, value=44000.0),
        ],
    )
    return model

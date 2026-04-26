"""EC 1.4.1.3.

Equilibrator
NADPH(aq) + NH3(aq) + 2-Oxoglutarate(aq) ⇌ H2O(l) + NADP(aq) + L-Glutamate(aq)
Keq = 7.2e5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_3s_2p
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


def add_glutamate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    nh4: str | None = None,
    oxoglutarate: str | None = None,
    glutamate: str | None = None,
    nadp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glutamate dehydrogenase (EC 1) to model."""
    rxn = default_name(rxn, n.glutamate_dehydrogenase)
    nadph = default_name(nadph, n.nadph)
    nh4 = default_name(nh4, n.nh4)
    oxoglutarate = default_name(oxoglutarate, n.oxoglutarate)
    glutamate = default_name(glutamate, n.glutamate)
    nadp = default_name(nadp, n.nadp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadph: -1.0,
                nh4: -1.0,
                oxoglutarate: -1.0,
                glutamate: 1.0,
                nadp: 1.0,
            },
        ),
        args=[
            nadph,
            nh4,
            oxoglutarate,
            glutamate,
            nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=104,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=1.54),
            default_kmp(model, rxn=rxn, par=kmp, value=0.64),
            default_keq(model, rxn=rxn, par=keq, value=7.2e5),
        ],
    )
    return model

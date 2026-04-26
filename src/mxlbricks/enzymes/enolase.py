"""EC 4.2.1.11.

Equilibrator
PEP(aq) + H2O(l) ⇌ 2-Phospho-D-glycerate(aq)
Keq = 0.2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_1s_1p
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


def add_enolase(
    model: Model,
    *,
    rxn: str | None = None,
    pep: str | None = None,
    pga2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add enolase (EC 4) to model."""
    rxn = default_name(rxn, n.enolase)
    pep = default_name(pep, n.pep)
    pga2 = default_name(pga2, n.pga2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                pep: -1.0,
                pga2: 1.0,
            },
        ),
        args=[
            pep,
            pga2,
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
            default_keq(model, rxn=rxn, par=keq, value=0.2),
        ],
    )
    return model

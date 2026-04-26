"""name.

EC 1.2.3.4

Equilibrator
Oxalate(aq) + O2(aq) + 2 H2O(l) ⇌ Hydrogen peroxide(aq) + 2 CO2(total)
Keq = 2.8e30 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_oxalate_oxidase(
    model: Model,
    *,
    rxn: str | None = None,
    oxalate: str | None = None,
    o2: str | None = None,
    h2o2: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add oxalate oxidase (name) to model."""
    rxn = default_name(rxn, n.oxalate_oxidase)
    oxalate = default_name(oxalate, n.oxalate)
    o2 = default_name(o2, n.o2)
    h2o2 = default_name(h2o2, n.h2o2)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                oxalate: -1.0,
                o2: -1.0,
                h2o2: 1.0,
                co2: 2.0,
            },
        ),
        args=[
            oxalate,
            o2,
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

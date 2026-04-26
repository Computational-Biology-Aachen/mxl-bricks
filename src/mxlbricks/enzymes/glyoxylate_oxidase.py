"""name.

EC 1.2.3.5

Equilibrator
Glyoxylate(aq) + H2O(l) + O2(aq) ⇌ Oxalate(aq) + H2O2(aq)
Keq = 2.5e28 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_glyoxylate_oxidase(
    model: Model,
    *,
    rxn: str | None = None,
    glyoxylate: str | None = None,
    o2: str | None = None,
    oxalate: str | None = None,
    h2o2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glyoxylate oxidase (name) to model."""
    rxn = default_name(rxn, n.glyoxylate_oxidase)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    o2 = default_name(o2, n.o2)
    oxalate = default_name(oxalate, n.oxalate)
    h2o2 = default_name(h2o2, n.h2o2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                # n.h2o(): -1.0,
                glyoxylate: -1.0,
                o2: -1.0,
                oxalate: 1.0,
                h2o2: 1.0,
            },
        ),
        args=[
            glyoxylate,
            o2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=1.0),
        ],
    )
    return model

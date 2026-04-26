"""catalase.

2 H2O2 <=> 2 H2O + O2

Equilibrator
2 H2O2(aq) ⇌ 2 H2O(l) + O2(aq)
Keq = 4.3e33 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_1s
from mxlbricks.utils import (
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_catalase(
    model: Model,
    *,
    rxn: str | None = None,
    h2o2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add catalase (catalase) to model."""
    rxn = default_name(rxn, n.catalase)
    h2o2 = default_name(h2o2, n.h2o2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry={
            h2o2: -1,
        },
        args=[
            h2o2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=760500.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=137.9),
        ],
    )
    return model

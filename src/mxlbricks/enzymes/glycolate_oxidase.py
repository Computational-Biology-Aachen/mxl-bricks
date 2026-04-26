"""glycolate oxidase.

O2 + Glycolate(chl) <=> H2O2 + Glyoxylate

Equilibrator
O2(aq) + Glycolate(aq) ⇌ H2O2(aq) + Glyoxylate(aq)
Keq = 3e15 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_1s, michaelis_menten_2s
from mxlbricks.utils import (
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_glycolate_oxidase_yokota(
    model: Model,
    *,
    rxn: str | None = None,
    glycolate: str | None = None,
    glyoxylate: str | None = None,
    h2o2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycolate oxidase without explicit oxygen concentration term."""
    rxn = default_name(rxn, n.glycolate_oxidase)
    glycolate = default_name(glycolate, n.glycolate)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    h2o2 = default_name(h2o2, n.h2o2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry={
            glycolate: -1,
            glyoxylate: 1,
            h2o2: 1,
        },
        args=[
            glycolate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=100,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.06),
        ],
    )
    return model


def add_glycolate_oxidase(
    model: Model,
    *,
    rxn: str | None = None,
    glycolate: str | None = None,
    glyoxylate: str | None = None,
    h2o2: str | None = None,
    o2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycolate oxidase (glycolate oxidase) to model."""
    rxn = default_name(rxn, n.glycolate_oxidase)
    glycolate = default_name(glycolate, n.glycolate)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    h2o2 = default_name(h2o2, n.h2o2)
    o2 = default_name(o2, n.o2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry={
            glycolate: -1,
            glyoxylate: 1,
            h2o2: 1,
        },
        args=[
            glycolate,
            o2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=100,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.06),
        ],
    )
    return model

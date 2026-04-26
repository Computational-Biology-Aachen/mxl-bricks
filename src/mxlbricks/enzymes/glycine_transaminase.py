"""glycine transaminase.

EC 2.6.1.4

Equilibrator
L-Glutamate(aq) + Glyoxylate(aq) ⇌ 2-Oxoglutarate(aq) + Glycine(aq)
Keq = 30 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    michaelis_menten_1s,
    michaelis_menten_2s,
    reversible_michaelis_menten_2s_2p,
)
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


def add_glycine_transaminase_yokota(
    model: Model,
    *,
    rxn: str | None = None,
    glyoxylate: str | None = None,
    glycine: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Yokota 1980 used reduced stoichiometry for the reaction."""
    rxn = default_name(rxn, n.glycine_transaminase)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    glycine = default_name(glycine, n.glycine)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glyoxylate: -1.0,
                glycine: 1.0,
            },
        ),
        args=[
            glyoxylate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=143.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=3.0),
        ],
    )
    return model


def add_glycine_transaminase_irreversible(
    model: Model,
    *,
    rxn: str | None = None,
    glutamate: str | None = None,
    glyoxylate: str | None = None,
    oxoglutarate: str | None = None,
    glycine: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycine transaminase irreversible (glycine transaminase) to model."""
    rxn = default_name(rxn, n.glycine_transaminase)
    glutamate = default_name(glutamate, n.glutamate)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    oxoglutarate = default_name(oxoglutarate, n.oxoglutarate)
    glycine = default_name(glycine, n.glycine)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glutamate: -1.0,
                glyoxylate: -1.0,
                oxoglutarate: 1.0,
                glycine: 1.0,
            },
        ),
        args=[
            glyoxylate,
            glutamate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=143.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=3.0),
        ],
    )

    return model


def add_glycine_transaminase(
    model: Model,
    *,
    rxn: str | None = None,
    glutamate: str | None = None,
    glyoxylate: str | None = None,
    oxoglutarate: str | None = None,
    glycine: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycine transaminase (glycine transaminase) to model."""
    rxn = default_name(rxn, n.glycine_transaminase)
    glutamate = default_name(glutamate, n.glutamate)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    oxoglutarate = default_name(oxoglutarate, n.oxoglutarate)
    glycine = default_name(glycine, n.glycine)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glutamate: -1.0,
                glyoxylate: -1.0,
                oxoglutarate: 1.0,
                glycine: 1.0,
            },
        ),
        args=[
            glyoxylate,
            glutamate,
            glycine,
            oxoglutarate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=143.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=3.0),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=30),
        ],
    )

    return model

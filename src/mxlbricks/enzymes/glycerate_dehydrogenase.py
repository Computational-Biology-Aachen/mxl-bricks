"""glycerate dehydrogenase.

NADH + Hydroxypyruvate <=> NAD  + D-Glycerate

Equilibrator
NADH(aq) + Hydroxypyruvate(aq) ⇌ NAD(aq) + D-Glycerate(aq)
Keq = 8.7e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    michaelis_menten_1s,
    reversible_michaelis_menten_2s_2p,
)
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_hpa_outflux(
    model: Model,
    *,
    rxn: str | None = None,
    hydroxypyruvate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add hpa outflux (glycerate dehydrogenase) to model."""
    rxn = default_name(rxn, n.glycerate_dehydrogenase)
    hydroxypyruvate = default_name(hydroxypyruvate, n.hydroxypyruvate)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry={
            hydroxypyruvate: -1.0,
        },
        args=[
            hydroxypyruvate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=398.0,  # Source
            ),
            static(model, n.kms(rxn), 0.12) if kms is None else kms,
        ],
    )

    return model


def add_glycerate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    hydroxypyruvate: str | None = None,
    nadh: str | None = None,
    glycerate: str | None = None,
    nad: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycerate dehydrogenase (glycerate dehydrogenase) to model."""
    rxn = default_name(rxn, n.glycerate_dehydrogenase)
    hydroxypyruvate = default_name(hydroxypyruvate, n.hydroxypyruvate)
    nadh = default_name(nadh, n.nadh)
    glycerate = default_name(glycerate, n.glycerate)
    nad = default_name(nad, n.nad)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadh: -1.0,
                hydroxypyruvate: -1.0,
                nad: 1.0,
                glycerate: 1.0,
            },
        ),
        args=[
            hydroxypyruvate,
            nadh,
            glycerate,
            nad,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=398.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.12),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=87000.0),
        ],
    )

    return model

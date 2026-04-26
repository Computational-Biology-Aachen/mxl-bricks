"""EC 1.1.99.14.

Equilibrator
------------
Glycolate(aq) + NAD (aq) ⇌ Glyoxylate(aq) + NADH(aq)
Keq = 6.2 x 10-8 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)

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


def add_glycolate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    glycolate: str | None = None,
    nad: str | None = None,
    glyoxylate: str | None = None,
    nadh: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycolate dehydrogenase (EC 1) to model."""
    rxn = default_name(rxn, n.glycolate_dehydrogenase)
    glycolate = default_name(glycolate, n.glycolate)
    nad = default_name(nad, n.nad)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    nadh = default_name(nadh, n.nadh)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycolate: -1,
                nad: -1,
                glyoxylate: 1,
                nadh: 1,
            },
        ),
        args=[
            glycolate,
            nad,
            glyoxylate,
            nadh,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.93,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.21),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=6.2e-8),
        ],
    )
    return model

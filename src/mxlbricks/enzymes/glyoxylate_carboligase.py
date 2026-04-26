"""glyoxylate carboligase == tartronate-semialdehyde synthase.

EC 4.1.1.47

Equilibrator
2 Glyoxylate + H2O <=> Tartronate semialdehyde + CO2(total)
Keq = 1.6e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_1s_2p
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


def add_glyoxylate_carboligase(
    model: Model,
    *,
    rxn: str | None = None,
    glyoxylate: str | None = None,
    tartronate_semialdehyde: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glyoxylate carboligase (glyoxylate carboligase == tartronate-semialdehyde synthase) to model."""
    rxn = default_name(rxn, n.glyoxylate_carboligase)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    tartronate_semialdehyde = default_name(
        tartronate_semialdehyde, n.tartronate_semialdehyde
    )
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glyoxylate: -2,
                tartronate_semialdehyde: 1,
                co2: 1,
            },
        ),
        args=[
            glyoxylate,
            tartronate_semialdehyde,
            co2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=18.9,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.9),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=1.6e4),
        ],
    )
    return model

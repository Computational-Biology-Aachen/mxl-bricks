"""name.

EC 1.17.1.9
    Metacyc: 1.2.1.2-RXN
    FORMATE + NAD ⇌ CARBON-DIOXIDE + NADH

    Equilibrator
    NAD (aq) + Formate(aq) + H2O(l) ⇌ NADH(aq) + CO2(total)
    Keq = 8.7e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_formate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    nad: str | None = None,
    formate: str | None = None,
    nadh: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add formate dehydrogenase (name) to model."""
    rxn = default_name(rxn, n.formate_dehydrogenase)
    nad = default_name(nad, n.nad)
    formate = default_name(formate, n.formate)
    nadh = default_name(nadh, n.nadh)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            stoichiometry={
                formate: -1.0,
                nad: -1.0,
                nadh: 1.0,
                co2: 1.0,
            },
        ),
        args=[
            nad,
            formate,
            nadh,
            co2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=2.9,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.011),
            default_kmp(model, rxn=rxn, par=kmp, value=0.18),
            default_keq(model, rxn=rxn, par=keq, value=8700.0),
        ],
    )
    return model

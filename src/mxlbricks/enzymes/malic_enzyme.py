"""malic enzyme == malate dehydrogenase decarboxylating.

EC 1.1.1.39

Equilibrator
    NAD (aq) + (S)-Malate(aq) + H2O(l) ⇌ NADH(aq) + Pyruvate(aq) + CO2(total)
    Keq = 0.2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_2s_3p
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


def add_malic_enzyme(
    model: Model,
    *,
    rxn: str | None = None,
    nad: str | None = None,
    malate: str | None = None,
    nadh: str | None = None,
    pyruvate: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add malic enzyme (malic enzyme == malate dehydrogenase decarboxylating) to model."""
    rxn = default_name(rxn, n.malic_enzyme)
    nad = default_name(nad, n.nad)
    malate = default_name(malate, n.malate)
    nadh = default_name(nadh, n.nadh)
    pyruvate = default_name(pyruvate, n.pyruvate)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nad: -1,
                malate: -1,
                nadh: 1,
                pyruvate: 1,
                co2: 1,
            },
        ),
        args=[
            nad,
            malate,
            nadh,
            pyruvate,
            co2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=39,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.003),
            default_kmp(model, rxn=rxn, par=kmp, value=0.00125),
            default_keq(model, rxn=rxn, par=keq, value=0.2),
        ],
    )
    return model

"""Tartronate semialdehyde(aq) + NADH(aq) ⇌ Glycerate(aq) + NAD (aq).

Keq = 1.6e5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25).
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


def add_tartronate_semialdehyde_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    tartronate_semialdehyde: str | None = None,
    nadh: str | None = None,
    glycerate: str | None = None,
    nad: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add tartronate semialdehyde reductase (Tartronate semialdehyde(aq) + NADH(aq) ⇌ Glycerate(aq) + NAD (aq)) to model."""
    rxn = default_name(rxn, n.tartronate_semialdehyde_reductase)
    tartronate_semialdehyde = default_name(
        tartronate_semialdehyde, n.tartronate_semialdehyde
    )
    nadh = default_name(nadh, n.nadh)
    glycerate = default_name(glycerate, n.glycerate)
    nad = default_name(nad, n.nad)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                tartronate_semialdehyde: -1,
                nadh: -1,
                glycerate: 1,
                nad: 1,
            },
        ),
        args=[
            tartronate_semialdehyde,
            nadh,
            glycerate,
            nad,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=243,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.05),
            default_kmp(model, rxn=rxn, par=kmp, value=0.28),
            default_keq(model, rxn=rxn, par=keq, value=1.6e5),
        ],
    )
    return model

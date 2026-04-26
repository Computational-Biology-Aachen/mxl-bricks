"""(2R,3S)-beta-Hydroxyasparate hydro-lyase (Iminosuccinate forming).

EC FIXME

Equilibrator
3-hydroxyaspartate(aq) ⇌ Iminoaspartate(aq) + H2O(l)
Keq = 4.0 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_1s_1p
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_hydroxyaspartate_hydrolase(
    model: Model,
    *,
    rxn: str | None = None,
    hydroxyaspartate: str | None = None,
    iminoaspartate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add hydroxyaspartate hydrolase ((2R,3S)-beta-Hydroxyasparate hydro-lyase (Iminosuccinate forming)) to model."""
    rxn = default_name(rxn, n.hydroxyaspartate_hydrolase)
    hydroxyaspartate = default_name(hydroxyaspartate, n.hydroxyaspartate)
    iminoaspartate = default_name(iminoaspartate, n.iminoaspartate)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry={
            hydroxyaspartate: -1,
            iminoaspartate: 1,
        },
        args=[
            hydroxyaspartate,
            iminoaspartate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.1),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=4.0),
        ],
    )
    return model

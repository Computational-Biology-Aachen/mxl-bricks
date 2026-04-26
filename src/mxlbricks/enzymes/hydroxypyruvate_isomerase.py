"""name.

EC 5.3.1.22

Equilibrator
Hydroxypyruvate(aq) ⇌ Tartronate semialdehyde(aq)
Keq = 0.5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_hydroxypyruvate_isomerase(
    model: Model,
    *,
    rxn: str | None = None,
    hydroxypyruvate: str | None = None,
    tartronate_semialdehyde: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add hydroxypyruvate isomerase (name) to model."""
    rxn = default_name(rxn, n.hydroxypyruvate_isomerase)

    hydroxypyruvate = default_name(hydroxypyruvate, n.hydroxypyruvate)
    tartronate_semialdehyde = default_name(
        tartronate_semialdehyde, n.tartronate_semialdehyde
    )

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                hydroxypyruvate: -1.0,
                tartronate_semialdehyde: 1.0,
            },
        ),
        args=[
            hydroxypyruvate,
            tartronate_semialdehyde,
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
            default_keq(model, rxn=rxn, par=keq, value=0.5),
        ],
    )
    return model

"""A5P <=> RU5P.

EC 5.3.1.13

Equilibrator
Arabinose-5-phosphate <=> Ru5P
Keq = 0.4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_a5p_isomerase(
    model: Model,
    *,
    rxn: str | None = None,
    a5p: str | None = None,
    ru5p: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add a5p isomerase (A5P <=> RU5P) to model."""
    rxn = default_name(rxn, n.a5p_isomerase)
    a5p = default_name(a5p, n.arabinose_5_phosphate)
    ru5p = default_name(ru5p, n.ru5p)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                a5p: -1.0,
                ru5p: 1.0,
            },
        ),
        args=[
            a5p,
            ru5p,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Clostridium tetani
                kcat_value=102,  # Clostridium tetani
            ),
            default_kms(model, rxn=rxn, par=kms, value=1.89),
            default_kmp(model, rxn=rxn, par=kmp, value=6.65),
            default_keq(model, rxn=rxn, par=keq, value=0.4),
        ],
    )
    return model

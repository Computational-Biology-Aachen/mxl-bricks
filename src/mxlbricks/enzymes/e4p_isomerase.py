"""ERYTHRULOSE_4P <=> ERYTHROSE_4P.

EC 5.3.1.34

Equilibrator
Erythrulose-4-phosphate(aq) ⇌ Erythrose-4-phosphate(aq)
Keq = 0.07 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_e4p_isomerase(
    model: Model,
    *,
    rxn: str | None = None,
    e4p: str | None = None,
    eu4p: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add e4p isomerase (ERYTHRULOSE_4P <=> ERYTHROSE_4P) to model."""
    rxn = default_name(rxn, n.e4p_isomerase)
    e4p = default_name(e4p, n.erythrulose_4p)
    eu4p = default_name(eu4p, n.e4p)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                e4p: -1.0,
                eu4p: 1.0,
            },
        ),
        args=[
            e4p,
            eu4p,
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
            default_keq(model, rxn=rxn, par=keq, value=0.07),
        ],
    )
    return model

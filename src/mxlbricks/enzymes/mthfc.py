"""Methenyltetrahydrofolate cyclohydrolase.

EC 3.5.4.9
Metacyc: MTHFC

10-FORMYL-THF_m + 0.07 PROTON_m <=> 5-10-METHENYL-THF_m + WATER_m

Equilibrator
10-Formyl-THF(aq) ⇌ 5,10-Methenyltetrahydrofolate(aq) + H2O(l)
Keq = 0.1 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_mthfc(
    model: Model,
    *,
    rxn: str | None = None,
    formyl_thf: str | None = None,
    methenyl_thf: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add mthfc (Methenyltetrahydrofolate cyclohydrolase) to model."""
    rxn = default_name(rxn, n.mthfc)
    formyl_thf = default_name(formyl_thf, n.formyl_thf)
    methenyl_thf = default_name(methenyl_thf, n.methenyl_thf)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                formyl_thf: -1,
                methenyl_thf: 1,
            },
        ),
        args=[
            formyl_thf,
            methenyl_thf,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=40.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.2),
            default_kmp(model, rxn=rxn, par=kmp, value=0.04),
            default_keq(model, rxn=rxn, par=keq, value=0.1),
        ],
    )

    return model

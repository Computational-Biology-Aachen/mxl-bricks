"""pyruvate dehydrogenase.

EC 1.2.4.1

Equilibrator
------------
    NAD (aq) + CoA(aq) + Pyruvate(aq) + H2O(l) ⇌ NADH(aq) + Acetyl-CoA(aq) + CO2(total)
    Keq = 2.6e7 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_3s_3p
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


def add_pyruvate_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    nad: str | None = None,
    coa: str | None = None,
    pyruvate: str | None = None,
    nadh: str | None = None,
    acetyl_coa: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add pyruvate dehydrogenase (pyruvate dehydrogenase) to model."""
    rxn = default_name(rxn, n.pyruvate_dehydrogenase)
    nad = default_name(nad, n.nad)
    coa = default_name(coa, n.coa)
    pyruvate = default_name(pyruvate, n.pyruvate)
    nadh = default_name(nadh, n.nadh)
    acetyl_coa = default_name(acetyl_coa, n.acetyl_coa)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nad: -1,
                coa: -1,
                pyruvate: -1,
                nadh: 1,
                acetyl_coa: 1,
                co2: 1,
            },
        ),
        args=[
            nad,
            coa,
            pyruvate,
            nadh,
            acetyl_coa,
            co2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=0.48,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.0124),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=2.6e7),
        ],
    )
    return model

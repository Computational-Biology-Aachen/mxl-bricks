"""Glycolyl-CoA + NADPH <=> Glycolaldehyde + NADP + CoA.

EC 1.2.1.12

Equilibrator

Keq = ? (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_glycolaldehyde_dehydrogenase(
    model: Model,
    *,
    rxn: str | None = None,
    glycolyl_coa: str | None = None,
    nadph: str | None = None,
    glycolaldehyde: str | None = None,
    nadp: str | None = None,
    coa: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycolaldehyde dehydrogenase (Glycolyl-CoA + NADPH <=> Glycolaldehyde + NADP + CoA) to model."""
    rxn = default_name(rxn, n.glycolaldehyde_dehydrogenase)
    glycolyl_coa = default_name(glycolyl_coa, n.glycolyl_coa)
    nadph = default_name(nadph, n.nadph)
    glycolaldehyde = default_name(glycolaldehyde, n.glycolaldehyde)
    nadp = default_name(nadp, n.nadp)
    coa = default_name(coa, n.coa)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycolyl_coa: -1.0,
                nadph: -1.0,
                glycolaldehyde: 1.0,
                nadp: 1.0,
                coa: 1.0,
            },
        ),
        args=[
            glycolyl_coa,
            nadph,
            glycolaldehyde,
            nadp,
            coa,
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
            default_keq(model, rxn=rxn, par=keq, value=1.0),
        ],
    )
    return model

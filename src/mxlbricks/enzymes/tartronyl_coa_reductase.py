"""TCR1: nadph + tarcoa -> nadp + coa + 2h3oppan.

Tartronyl-Coa + NADPH -> Tartronate-semialdehyde + NADP + CoA
dG' = 29.78
Keq = 6.06e-6
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


def add_tartronyl_coa_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    tartronyl_coa: str | None = None,
    nadp: str | None = None,
    tartronate_semialdehyde: str | None = None,
    coa: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add tartronyl coa reductase (TCR1: nadph + tarcoa -> nadp + coa + 2h3oppan) to model."""
    rxn = default_name(rxn, n.tartronyl_coa_reductase)
    nadph = default_name(nadph, n.nadph)
    tartronyl_coa = default_name(tartronyl_coa, n.tartronyl_coa)
    nadp = default_name(nadp, n.nadp)
    tartronate_semialdehyde = default_name(
        tartronate_semialdehyde, n.tartronate_semialdehyde
    )
    coa = default_name(coa, n.coa)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadph: -1,
                tartronyl_coa: -1,
                nadp: 1,
                tartronate_semialdehyde: 1,
                coa: 1,
            },
        ),
        args=[
            # substrates
            tartronyl_coa,
            nadph,
            # products
            tartronate_semialdehyde,
            nadp,
            coa,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.4,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.03),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=6.06e-06),
        ],
    )
    return model

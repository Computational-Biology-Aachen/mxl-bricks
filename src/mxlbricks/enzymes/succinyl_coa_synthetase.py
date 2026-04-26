"""name.

EC 6.2.1.5

Metacyc: SUCCCOASYN-RXN
SUC_m + CO-A_m + ATP_m <=> SUC-COA_m + ADP_m + Pi_m

Equilibrator
Succinate(aq) + CoA(aq) + ATP(aq) ⇌ Succinyl-CoA(aq) + ADP(aq) + Pi(aq)
Keq = 2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_succinyl_coa_synthetase(
    model: Model,
    *,
    rxn: str | None = None,
    succinate: str | None = None,
    coa: str | None = None,
    succinyl_coa: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add succinyl coa synthetase (name) to model."""
    rxn = default_name(rxn, n.succinyl_coa_synthetase)
    succinate = default_name(succinate, n.succinate)
    coa = default_name(coa, n.coa)
    succinyl_coa = default_name(succinyl_coa, n.succinyl_coa)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                succinate: -1,
                coa: -1,
                succinyl_coa: 1,
                adp: 1,
                pi: 1,
            },
        ),
        args=[
            succinate,
            coa,
            succinyl_coa,
            adp,
            pi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=44.73,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.25),
            default_kmp(model, rxn=rxn, par=kmp, value=0.041),
            default_keq(model, rxn=rxn, par=keq, value=2.0),
        ],
    )

    return model

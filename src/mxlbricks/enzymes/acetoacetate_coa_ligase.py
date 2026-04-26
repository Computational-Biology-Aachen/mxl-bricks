"""acetoacetate_coa_ligase.

EC 6.2.1.16

Metacyc (ACETOACETATE--COA-LIGASE-RXN):
    3-KETOBUTYRATE_m + ATP_m + CO-A_m
    --> ACETOACETYL-COA_m + AMP_m + Diphosphate_m + 0.92 PROTON_m

Equilibrator
    Acetoacetate(aq) + ATP(aq) + CoA(aq) ⇌ Acetoacetyl-CoA(aq) + AMP(aq) + Diphosphate(aq)
    Keq = 2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_acetoacetate_coa_ligase(
    model: Model,
    *,
    rxn: str | None = None,
    acac: str | None = None,
    atp: str | None = None,
    coa: str | None = None,
    acac_coa: str | None = None,
    amp: str | None = None,
    ppi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add acetoacetate coa ligase (acetoacetate_coa_ligase) to model."""
    rxn = default_name(rxn, n.acetoacetate_coa_ligase)
    acac = default_name(acac, n.acetoacetate)
    atp = default_name(atp, n.atp)
    coa = default_name(coa, n.coa)
    acac_coa = default_name(acac_coa, n.acetoacetyl_coa)
    amp = default_name(amp, n.amp)
    ppi = default_name(ppi, n.ppi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                acac: -1,
                atp: -1,
                coa: -1,
                acac_coa: 1,
                amp: 1,
                ppi: 1,
            },
        ),
        args=[
            acac,
            atp,
            coa,
            acac_coa,
            amp,
            ppi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=5.89,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.07),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=2.0),
        ],
    )

    return model

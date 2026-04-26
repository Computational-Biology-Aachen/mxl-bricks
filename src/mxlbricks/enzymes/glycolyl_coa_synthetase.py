"""GCS: atp + coa + glyclt -> Diphosphate + amp + glyccoa.

dG' = 9.25
Keq = 0.024
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    michaelis_menten_3s,
    reversible_michaelis_menten_3s_3p,
)
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


def add_glycolyl_coa_synthetase_irrev(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    coa: str | None = None,
    glycolate: str | None = None,
    glycolyl_coa: str | None = None,
    ppi: str | None = None,
    amp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycolyl coa synthetase irrev (GCS: atp + coa + glyclt -> Diphosphate + amp + glyccoa) to model."""
    rxn = default_name(rxn, n.glycolyl_coa_synthetase)
    atp = default_name(atp, n.atp)
    coa = default_name(coa, n.coa)
    glycolate = default_name(glycolate, n.glycolate)
    glycolyl_coa = default_name(glycolyl_coa, n.glycolyl_coa)
    ppi = default_name(ppi, n.ppi)
    amp = default_name(amp, n.amp)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_3s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                atp: -1,
                coa: -1,
                glycolate: -1,
                glycolyl_coa: 1,
                ppi: 1,
                amp: 1,
            },
        ),
        args=[
            atp,
            coa,
            glycolate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=4.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=13.0),
        ],
    )
    return model


def add_glycolyl_coa_synthetase(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    coa: str | None = None,
    glycolate: str | None = None,
    glycolyl_coa: str | None = None,
    ppi: str | None = None,
    amp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycolyl coa synthetase (GCS: atp + coa + glyclt -> Diphosphate + amp + glyccoa) to model."""
    rxn = default_name(rxn, n.glycolyl_coa_synthetase)
    atp = default_name(atp, n.atp)
    coa = default_name(coa, n.coa)
    glycolate = default_name(glycolate, n.glycolate)
    glycolyl_coa = default_name(glycolyl_coa, n.glycolyl_coa)
    ppi = default_name(ppi, n.ppi)
    amp = default_name(amp, n.amp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                atp: -1,
                coa: -1,
                glycolate: -1,
                glycolyl_coa: 1,
                ppi: 1,
                amp: 1,
            },
        ),
        args=[
            atp,
            coa,
            glycolate,
            glycolyl_coa,
            ppi,
            amp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=4.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=13.0),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=0.024),
        ],
    )
    return model

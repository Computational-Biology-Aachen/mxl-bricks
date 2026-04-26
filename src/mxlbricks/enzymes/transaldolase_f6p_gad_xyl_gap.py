"""F6P + Glycolaldehyde <=> GAP + XYLULOSE.

EC 2.2.1.2

Equilibrator
Fructose-6-phosphate(aq) + Glycolaldehyde(aq)
    ⇌ Glyceraldehyde 3-phosphate(aq) + Xylulose(aq)
Keq = 4.8e-4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_2s_2p
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


def add_transaldolase_f6p_gad_xyl_gap(
    model: Model,
    *,
    rxn: str | None = None,
    f6p: str | None = None,
    glycolaldehyde: str | None = None,
    gap: str | None = None,
    xylulose: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add transaldolase f6p gad xyl gap (F6P + Glycolaldehyde <=> GAP + XYLULOSE) to model."""
    rxn = default_name(rxn, n.transaldolase_f6p_gad_gap_xyl)
    f6p = default_name(f6p, n.f6p)
    glycolaldehyde = default_name(glycolaldehyde, n.glycolaldehyde)
    gap = default_name(gap, n.gap)
    xylulose = default_name(xylulose, n.xylulose)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                f6p: -1.0,
                glycolaldehyde: -1.0,
                gap: 1.0,
                xylulose: 1.0,
            },
        ),
        args=[
            f6p,
            glycolaldehyde,
            gap,
            xylulose,
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
            default_keq(model, rxn=rxn, par=keq, value=4.8e-4),
        ],
    )
    return model

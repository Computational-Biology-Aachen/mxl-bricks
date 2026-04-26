"""Lumped reaction of Glyceraldehyde 3-phosphate dehydrogenase (GADPH) and Phosphoglycerate kinase (PGK).

    3-Phospho-D-glycerate(aq) + ATP(aq) ⇌ 3-Phospho-D-glyceroyl phosphate(aq) + ADP(aq)
    3-Phospho-D-glyceroyl phosphate(aq) + NADPH(aq) ⇌ D-Glyceraldehyde 3-phosphate(aq) + NADP (aq) + Orthophosphate(aq)
Into
    3-Phospho-D-glycerate(aq) + ATP(aq) + NADPH(aq) ⇌ D-Glyceraldehyde 3-phosphate(aq) + ADP(aq) + Orthophosphate(aq) + NADP(aq).

Equilibrator
    Keq = 6.0e-4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_3s_4p
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


def lumped_pgk_gadph(
    model: Model,
    *,
    rxn: str | None = None,
    pga: str | None = None,
    atp: str | None = None,
    nadph: str | None = None,
    gap: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    nadp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add lumped pgk gadph (Lumped reaction of Glyceraldehyde 3-phosphate dehydrogenase (GADPH) and Phosphoglycerate kinase (PGK)) to model."""
    rxn = default_name(rxn, n.pgk_gadph)
    pga = default_name(pga, n.pga)
    atp = default_name(atp, n.atp)
    nadph = default_name(nadph, n.nadph)
    gap = default_name(gap, n.gap)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)
    nadp = default_name(nadp, n.nadp)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_4p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                pga: -1.0,
                atp: -1.0,
                nadph: -1.0,
                gap: 1.0,
                adp: 1.0,
                pi: 1.0,
                nadp: 1.0,
            },
        ),
        args=[
            pga,
            atp,
            nadph,
            gap,
            adp,
            pi,
            nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=537,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.18),
            default_kmp(model, rxn=rxn, par=kmp, value=0.27),
            default_keq(model, rxn=rxn, par=keq, value=6.0e-4),
        ],
    )
    return model

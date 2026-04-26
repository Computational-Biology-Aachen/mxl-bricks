"""Glycolaldehyde + S7P <=> RIBOSE_5P  + ERYTHRULOSE.

EC 2.2.1.1

Equilibrator
Glycolaldehyde(aq) + Sedoheptulose-7-phosphate(aq)
    ⇌ Ribose-5-phosphate(aq) + Erythrulose(aq)
Keq = 0.5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
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


def add_transketolase_gad_s7p_eru_r5p(
    model: Model,
    *,
    rxn: str | None = None,
    glycolaldehyde: str | None = None,
    s7p: str | None = None,
    r5p: str | None = None,
    erythrulose: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add transketolase gad s7p eru r5p (Glycolaldehyde + S7P <=> RIBOSE_5P  + ERYTHRULOSE) to model."""
    rxn = default_name(rxn, n.transketolase_gad_s7p_r5p_eru)
    glycolaldehyde = default_name(glycolaldehyde, n.glycolaldehyde)
    s7p = default_name(s7p, n.s7p)
    r5p = default_name(r5p, n.r5p)
    erythrulose = default_name(erythrulose, n.erythrulose)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycolaldehyde: -1.0,
                s7p: -1.0,
                r5p: 1.0,
                erythrulose: 1.0,
            },
        ),
        args=[
            glycolaldehyde,
            s7p,
            r5p,
            erythrulose,
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
            default_keq(model, rxn=rxn, par=keq, value=0.5),
        ],
    )
    return model

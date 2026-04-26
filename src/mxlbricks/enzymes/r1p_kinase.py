"""R1P + ATP  <=> RUBP + ADP.

EC FIXME

Equilibrator
Ribose-1-phosphate(aq) + ATP(aq) ⇌ Ribulose-1,5-bisphosphate(aq) + ADP(aq)
Keq = 4.4e6 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_2s
from mxlbricks.utils import (
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_r1p_kinase(
    model: Model,
    *,
    rxn: str | None = None,
    r1p: str | None = None,
    atp: str | None = None,
    rubp: str | None = None,
    adp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add r1p kinase (R1P + ATP  <=> RUBP + ADP) to model."""
    rxn = default_name(rxn, n.r1p_kinase)
    r1p = default_name(r1p, n.r1p)
    atp = default_name(atp, n.atp)
    rubp = default_name(rubp, n.rubp)
    adp = default_name(adp, n.adp)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                r1p: -1.0,
                atp: -1.0,
                rubp: 1.0,
                adp: 1.0,
            },
        ),
        args=[
            r1p,
            atp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.1),
        ],
    )
    return model

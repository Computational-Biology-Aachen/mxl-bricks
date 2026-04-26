"""SBPase.

EC 3.1.3.37

Equilibrator
    H2O(l) + Sedoheptulose 1,7-bisphosphate(aq)
    ⇌ Orthophosphate(aq) + Sedoheptulose 7-phosphate(aq)
    Keq = 2e2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_1s_1i
from mxlbricks.utils import (
    default_kis,
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_sbpase(
    model: Model,
    *,
    rxn: str | None = None,
    sbp: str | None = None,
    s7p: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    ki: str | None = None,
) -> Model:
    """Add sbpase (SBPase) to model."""
    rxn = default_name(rxn, n.sbpase)
    sbp = default_name(sbp, n.sbp)
    s7p = default_name(s7p, n.s7p)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s_1i,
        stoichiometry={
            sbp: -1,
            s7p: 1,
        },
        args=[
            sbp,
            pi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=0.04 * 8,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.013),
            default_kis(model, rxn=rxn, par=ki, substrate=n.pi(), value=12.0),
        ],
    )
    return model

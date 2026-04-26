"""fructose-1,6-bisphosphatase.

EC 3.1.3.11

Equilibrator

Equilibrator
    H2O(l) + D-Fructose 1,6-bisphosphate(aq) ⇌ Orthophosphate(aq) + D-Fructose 6-phosphate(aq)
    Keq = 1.2e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import michaelis_menten_1s_2i
from mxlbricks.utils import (
    default_kis,
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_fbpase(
    model: Model,
    *,
    rxn: str | None = None,
    fbp: str | None = None,
    f6p: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    ki_f6p: str | None = None,
    ki_pi: str | None = None,
) -> Model:
    """Add fbpase (fructose-1,6-bisphosphatase) to model."""
    rxn = default_name(rxn, n.fbpase)
    fbp = default_name(fbp, n.fbp)
    f6p = default_name(f6p, n.f6p)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s_2i,
        stoichiometry={
            fbp: -1,
            f6p: 1,
        },
        args=[
            fbp,
            f6p,
            pi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=0.2 * 8,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.03),
            default_kis(model, par=ki_f6p, rxn=rxn, substrate=f6p, value=0.7),
            default_kis(model, par=ki_pi, rxn=rxn, substrate=pi, value=12.0),
        ],
    )
    return model

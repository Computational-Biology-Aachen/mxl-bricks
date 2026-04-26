"""Zeaxanthin Epoxidase (stroma):.

Zeaxanthin + NADPH + O2 -> Anteraxanthin + NADP + H2O
Antheraxanthin + NADPH + O2 -> Violaxanthin + NADP + H2O.
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_name,
    default_par,
)


def add_zeaxanthin_epoxidase(
    model: Model,
    *,
    rxn: str | None = None,
    vx: str | None = None,
    zx: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add zeaxanthin epoxidase (Zeaxanthin Epoxidase (stroma):) to model."""
    rxn = default_name(rxn, n.zeaxanthin_epoxidase)
    vx = default_name(vx, n.vx)
    zx = default_name(zx, n.zx)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry={
            vx: 1,
        },
        args=[
            zx,
            default_par(model, par=kf, name=n.kf(rxn), value=0.00024),
        ],
    )
    return model

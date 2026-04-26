"""Enzyme building block: ex atp."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_kf,
    default_name,
)


def add_atp_consumption(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add atp consumption to model."""
    rxn = default_name(rxn, n.ex_atp)
    atp = default_name(atp, n.atp)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry={
            atp: -1,
        },
        args=[
            atp,
            default_kf(model, rxn=rxn, par=kf, value=1.0),
        ],
    )
    return model

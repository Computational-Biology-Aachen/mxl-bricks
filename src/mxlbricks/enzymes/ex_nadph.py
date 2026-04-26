"""name.

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_kf,
    default_name,
)


def add_nadph_consumption(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add nadph consumption (name) to model."""
    rxn = default_name(rxn, n.ex_nadph)
    nadph = default_name(nadph, n.nadph)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry={
            nadph: -1,
        },
        args=[
            nadph,
            default_kf(model, rxn=rxn, par=kf, value=1.0),
        ],
    )
    return model

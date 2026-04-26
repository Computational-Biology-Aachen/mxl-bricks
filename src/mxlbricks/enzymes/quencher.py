"""Enzyme building block: quencher."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_name,
    filter_stoichiometry,
    static,
)


def add_quenching_reaction(
    model: Model,
    *,
    rxn: str | None = None,
    energy: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add quenching reaction to model."""
    rxn = default_name(rxn, n.quencher)
    energy = default_name(energy, n.energy)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                energy: -1.0,
            },
        ),
        args=[
            energy,
            static(model, kf := n.kre(rxn), 1.0) if kf is None else kf,
        ],
    )
    return model

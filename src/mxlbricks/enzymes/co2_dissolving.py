"""Enzyme building block: co2 dissolving."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import diffusion
from mxlbricks.utils import (
    default_kf,
    default_name,
)


def add_co2_dissolving(
    model: Model,
    *,
    rxn: str | None = None,
    co2: str | None = None,
    co2_atmosphere: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add co2 dissolving to model."""
    rxn = default_name(rxn, n.co2_dissolving)
    co2 = default_name(co2, n.co2)
    co2_atmosphere = default_name(co2_atmosphere, n.co2_atmosphere)

    model.add_reaction(
        name=rxn,
        fn=diffusion,
        stoichiometry={
            co2: 1,
        },
        args=[
            co2,
            co2_atmosphere,
            default_kf(model, rxn=rxn, par=kf, value=4.5),
        ],
    )
    return model

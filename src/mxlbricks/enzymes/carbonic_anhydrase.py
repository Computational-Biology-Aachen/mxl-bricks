"""name.

EC 4.2.1.1

Equilibrator

hco3:co2 is ~50:1 according to Straßburger
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import reversible_mass_action_keq_1s_1p
from mxlbricks.utils import (
    default_keq,
    default_kf,
    default_name,
)


def add_carbonic_anhydrase_mass_action(
    model: Model,
    *,
    rxn: str | None = None,
    s1: str | None = None,
    p1: str | None = None,
    kf: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add carbonic anhydrase mass action (name) to model."""
    rxn = default_name(rxn, n.carbonic_anhydrase)
    s1 = default_name(s1, n.co2)
    p1 = default_name(p1, n.hco3)

    model.add_reaction(
        name=rxn,
        fn=reversible_mass_action_keq_1s_1p,
        stoichiometry={
            s1: -1,
            p1: 1,
        },
        args=[
            s1,
            p1,
            default_kf(model, rxn=rxn, par=kf, value=1000),
            default_keq(model, rxn=rxn, par=keq, value=50),
        ],
    )
    return model

"""Plastid terminal oxidase.

2 QH2 + O2 -> 2 Q + 2 H2O
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_2s
from mxlbricks.utils import (
    default_name,
    static,
)


def add_ptox(
    model: Model,
    *,
    rxn: str | None = None,
    pq_ox: str | None = None,
    pq_red: str | None = None,
    o2: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add ptox (Plastid terminal oxidase) to model."""
    rxn = default_name(rxn, n.ptox)
    pq_ox = default_name(pq_ox, n.pq_ox)
    pq_red = default_name(pq_red, n.pq_red)
    o2 = default_name(o2, lambda: n.o2("_lumen"))

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry={
            pq_ox: 1,
        },
        args=[
            pq_red,
            o2,
            static(model, "kPTOX", 0.01) if kf is None else kf,
        ],
    )
    return model

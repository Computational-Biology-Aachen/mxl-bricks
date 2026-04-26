"""NAD(P)H dehydrogenase-like complex (NDH).

PQH2 -> PQ

"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_kf,
    default_name,
)


def add_ndh(
    model: Model,
    *,
    rxn: str | None = None,
    pq_ox: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add ndh (NAD(P)H dehydrogenase-like complex (NDH)) to model."""
    rxn = default_name(rxn, n.ndh)
    pq_ox = default_name(pq_ox, n.pq_ox)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry={
            pq_ox: -1,
        },
        args=[
            pq_ox,
            default_kf(model, par=kf, rxn=rxn, value=0.002),
        ],
    )
    return model

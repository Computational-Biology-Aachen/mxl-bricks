"""name.

EC FIXME

Equilibrator
Monodehydroascorbate(aq) + 0.5 NAD (aq) ⇌ Dehydroascorbate(aq) + 0.5 NADH(aq)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_kf,
    default_name,
)


def _rate_mda_reductase(
    mda: float,
    k3: float,
) -> float:
    return k3 * mda**2


def add_mda_reductase1(
    model: Model,
    *,
    rxn: str | None = None,
    mda: str | None = None,
    dha: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add mda reductase1 (name) to model."""
    rxn = default_name(rxn, n.mda_reductase1)
    mda = default_name(mda, n.mda)
    dha = default_name(dha, n.dha)

    model.add_reaction(
        name=rxn,
        fn=_rate_mda_reductase,
        stoichiometry={
            mda: -2,
            dha: 1,
        },
        args=[
            mda,
            default_kf(model, rxn=rxn, par=kf, value=500.0),
        ],
    )
    return model

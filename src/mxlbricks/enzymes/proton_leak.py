"""name.

EC FIXME

Equilibrator
"""

from mxlpy import Derived, Model

from mxlbricks import names as n
from mxlbricks.fns import protons_stroma_2016
from mxlbricks.utils import (
    default_name,
    static,
)


def _neg_one_div_by(x: float) -> float:
    """Return -1/x; used for negated unit stoichiometry scaled by buffering capacity."""
    return -1.0 / x


def _rate_leak(
    protons_lumen: float,
    ph_stroma: float,
    k_leak: float,
) -> float:
    """Passive proton leak across the thylakoid membrane, proportional to the proton gradient."""
    return k_leak * (protons_lumen - protons_stroma_2016(ph_stroma))


def add_proton_leak(
    model: Model,
    *,
    rxn: str | None = None,
    h_lumen: str | None = None,
    ph_stroma: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add proton leak (name) to model."""
    rxn = default_name(rxn, n.proton_leak)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    ph_stroma = default_name(ph_stroma, n.ph)

    model.add_reaction(
        name=rxn,
        fn=_rate_leak,
        stoichiometry={
            h_lumen: Derived(fn=_neg_one_div_by, args=["bH"]),
        },
        args=[
            h_lumen,
            ph_stroma,
            static(model, n.kf(rxn), 10.0) if kf is None else kf,
        ],
    )
    return model

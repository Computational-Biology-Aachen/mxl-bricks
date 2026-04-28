"""name.

EC FIXME

Equilibrator
        2-Oxoglutarate(aq) + ATP(aq) + 2 ferredoxin(red)(aq) + NH4 (aq)
        ⇌ Glutamate(aq) + ADP(aq) + 2 ferredoxin(ox)(aq) + Orthophosphate(aq)
    K'eq = 2.4e13
"""

from mxlpy import Derived, Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_kf,
    default_name,
    filter_stoichiometry,
    static,
)


def _two_times_convf(convf: float) -> float:
    """Return 2 * convf for scaling Fd_ox stoichiometry."""
    return 2.0 * convf


def _rate_nitrogen_fixation(
    oxo: float,
    atp: float,
    fd_red: float,
    nh4: float,
    k_fwd: float,
    convf: float,
) -> float:
    """Nitrogen fixation rate: mass action on 2-oxoglutarate, ATP, NH4+ and reduced ferredoxin."""
    return k_fwd * oxo * atp * nh4 * (2 * fd_red * convf)


def add_nitrogen_metabolism(
    model: Model,
    *,
    rxn: str | None = None,
    oxoglutarate: str | None = None,
    atp: str | None = None,
    fd_red: str | None = None,
    nh4: str | None = None,
    glutamate: str | None = None,
    fd_ox: str | None = None,
    kf: str | None = None,
    convf: str | None = None,
) -> Model:
    """Add nitrogen metabolism (name) to model."""
    rxn = default_name(rxn, n.nitrogen_fixation)
    oxoglutarate = default_name(oxoglutarate, n.oxoglutarate)
    atp = default_name(atp, n.atp)
    fd_red = default_name(fd_red, n.fd_red)
    nh4 = default_name(nh4, n.nh4)
    glutamate = default_name(glutamate, n.glutamate)
    fd_ox = default_name(fd_ox, n.fd_ox)

    convf = static(model, n.convf(), 3.2e-2) if convf is None else convf

    model.add_reaction(
        rxn,
        _rate_nitrogen_fixation,
        stoichiometry=filter_stoichiometry(
            model,
            {
                atp: -1.0,  # mM
                nh4: -1.0,  # mM
                glutamate: 1.0,  # mM
                fd_ox: Derived(fn=_two_times_convf, args=[convf]),
            },
        ),
        args=[
            oxoglutarate,
            atp,
            fd_red,
            nh4,
            default_kf(model, par=kf, rxn=rxn, value=1.0),
            convf,
        ],
    )

    return model

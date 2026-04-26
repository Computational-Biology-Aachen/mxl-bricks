"""Enzyme building block: fd reductase."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    default_vmax,
    filter_stoichiometry,
)


def _rate_ferredoxin_reductase(
    Fd: float,
    Fdred: float,
    A1: float,
    A2: float,
    kFdred: float,
    Keq_FAFd: float,
) -> float:
    """Rate of the redcution of Fd by the activity of PSI.

    used to be equall to the rate of PSI but now
    alternative electron pathway from Fd allows for the production of ROS
    hence this rate has to be separate.
    """
    return kFdred * Fd * A1 - kFdred / Keq_FAFd * Fdred * A2


def add_ferredoxin_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    fd_ox: str | None = None,
    fd_red: str | None = None,
    a1: str | None = None,
    a2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    keq: str,  # derived from PSI
) -> Model:
    """Add ferredoxin reductase to model."""
    rxn = default_name(rxn, n.ferredoxin_reductase)
    fd_ox = default_name(fd_ox, n.fd_ox)
    fd_red = default_name(fd_red, n.fd_red)
    a1 = default_name(a1, n.a1)
    a2 = default_name(a2, n.a2)

    model.add_reaction(
        name=rxn,
        fn=_rate_ferredoxin_reductase,
        stoichiometry=filter_stoichiometry(
            model,
            {
                fd_ox: -1,
                fd_red: 1,
            },
        ),
        args=[
            fd_ox,
            fd_red,
            a1,
            a2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=2.5e5,  # Source
            ),
            keq,  # no default value, derived from PSI
        ],
    )
    return model

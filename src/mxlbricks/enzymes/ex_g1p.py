"""name.

Equilibrator
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    default_vmax,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def _rate_starch(
    g1p: float,
    atp: float,
    adp: float,
    pi: float,
    pga: float,
    f6p: float,
    fbp: float,
    v_st: float,
    kmst1: float,
    kmst2: float,
    ki_st: float,
    kast1: float,
    kast2: float,
    kast3: float,
) -> float:
    """Starch synthesis rate via G1P+ATP with ADP inhibition and allosteric activation by PGA/F6P/FBP."""
    return (
        v_st
        * g1p
        * atp
        / (
            (g1p + kmst1)
            * (
                (1 + adp / ki_st) * (atp + kmst2)
                + kmst2 * pi / (kast1 * pga + kast2 * f6p + kast3 * fbp)
            )
        )
    )


def add_g1p_efflux(
    model: Model,
    *,
    rxn: str | None = None,
    g1p: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    pga: str | None = None,
    f6p: str | None = None,
    fbp: str | None = None,
    starch: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_g1p: str | None = None,
    km_atp: str | None = None,
    ki: str | None = None,
    ka_pga: str | None = None,
    ka_f6p: str | None = None,
    ka_fbp: str | None = None,
) -> Model:
    """Add g1p efflux (name) to model."""
    rxn = default_name(rxn, n.ex_g1p)
    g1p = default_name(g1p, n.g1p)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)
    pga = default_name(pga, n.pga)
    f6p = default_name(f6p, n.f6p)
    fbp = default_name(fbp, n.fbp)
    starch = default_name(starch, n.starch)

    km_g1p = static(model, n.km(rxn, n.g1p()), 0.08) if km_g1p is None else km_g1p
    km_atp = static(model, n.km(rxn, n.atp()), 0.08) if km_atp is None else km_atp
    ki = static(model, n.ki(rxn), 10.0) if ki is None else ki
    ka_pga = static(model, n.ka(rxn, n.pga()), 0.1) if ka_pga is None else ka_pga
    ka_f6p = static(model, n.ka(rxn, n.f6p()), 0.02) if ka_f6p is None else ka_f6p
    ka_fbp = static(model, n.ka(rxn, n.fbp()), 0.02) if ka_fbp is None else ka_fbp

    model.add_reaction(
        name=rxn,
        fn=_rate_starch,
        stoichiometry=filter_stoichiometry(
            model,
            {
                n.g1p(): -1.0,
                n.atp(): -1.0,
                n.adp(): 1.0,
            },
            optional={
                n.starch(): 1.0,
            },
        ),
        args=[
            n.g1p(),
            n.atp(),
            n.adp(),
            n.pi(),
            n.pga(),
            n.f6p(),
            n.fbp(),
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=0.04 * 8,  # Source
            ),
            km_g1p,
            km_atp,
            ki,
            ka_pga,
            ka_f6p,
            ka_fbp,
        ],
    )
    return model

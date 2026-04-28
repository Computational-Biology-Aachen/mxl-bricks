"""Ferredoxin-NADP reductase.

2 reduced ferredoxin + NADP+ + H+ ⇌ 2 oxidized ferredoxin + NADPH

EC 1.18.1.2

Equilibrator
"""

import math

from mxlpy import Derived, Model

from mxlbricks import names as n
from mxlbricks.fns import (
    michaelis_menten_1s,
    michaelis_menten_2s,
    value,
)
from mxlbricks.utils import (
    default_km,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
    static,
)


def _keq_fnr(
    E0_Fd: float,
    F: float,
    E0_NADP: float,
    pHstroma: float,
    dG_pH: float,
    RT: float,
) -> float:
    """Equilibrium constant for FNR: Fd-mediated NADP+ reduction, pH-corrected."""
    dg1 = -E0_Fd * F
    dg2 = -2 * E0_NADP * F
    dg = -2 * dg1 + dg2 + dG_pH * pHstroma
    return math.exp(-dg / RT)


def _rate_fnr2016(
    fd_ox: float,
    fd_red: float,
    nadph: float,
    nadp: float,
    vmax: float,
    km_fd_red: float,
    km_nadph: float,
    keq: float,
) -> float:
    """FNR rate (2016 formulation): reversible ping-pong with Fd^2 stoichiometry, mmol/mmol_Chl units."""
    fdred = fd_red / km_fd_red
    fdox = fd_ox / km_fd_red
    nadph = nadph / km_nadph
    nadp = nadp / km_nadph
    return (
        vmax
        * (fdred**2 * nadp - fdox**2 * nadph / keq)
        / ((1 + fdred + fdred**2) * (1 + nadp) + (1 + fdox + fdox**2) * (1 + nadph) - 1)
    )


def _rate_fnr_2019(
    Fd_ox: float,
    Fd_red: float,
    NADPH: float,
    NADP: float,
    KM_FNR_F: float,
    KM_FNR_N: float,
    vmax: float,
    Keq_FNR: float,
    convf: float,
) -> float:
    """FNR rate (2019 formulation): same as 2016 but NADP/H concentrations scaled by convf."""
    fdred = Fd_red / KM_FNR_F
    fdox = Fd_ox / KM_FNR_F
    nadph = NADPH / convf / KM_FNR_N
    nadp = NADP / convf / KM_FNR_N
    return (
        vmax
        * (fdred**2 * nadp - fdox**2 * nadph / Keq_FNR)
        / ((1 + fdred + fdred**2) * (1 + nadp) + (1 + fdox + fdox**2) * (1 + nadph) - 1)
    )


def add_fnr_mmol_chl(
    model: Model,
    *,
    rxn: str | None = None,
    fd_ox: str | None = None,
    fd_red: str | None = None,
    nadph: str | None = None,
    nadp: str | None = None,
    ph: str | None = None,
    kcat: str | None = None,
    km_fd_red: str | None = None,
    km_nadp: str | None = None,
    e0: str | None = None,
) -> Model:
    """Add fnr mmol chl (Ferredoxin-NADP reductase) to model."""
    rxn = default_name(rxn, n.fnr)
    fd_ox = default_name(fd_ox, n.fd_ox)
    fd_red = default_name(fd_red, n.fd_red)
    nadph = default_name(nadph, n.nadph)
    nadp = default_name(nadp, n.nadp)
    ph = default_name(ph, n.ph)

    model.add_derived(
        keq := n.keq(rxn),
        fn=_keq_fnr,
        args=[
            "E^0_Fd",
            "F",
            "E^0_NADP",
            ph,
            "dG_pH",
            "RT",
        ],
    )

    model.add_reaction(
        name=rxn,
        fn=_rate_fnr2016,
        stoichiometry={
            fd_ox: 2,
        },
        args=[
            fd_ox,
            fd_red,
            nadph,
            nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=3.0,  # Source
                kcat_value=500.0,  # Source
            ),
            default_km(model, par=km_fd_red, rxn=rxn, subs=fd_red, value=1.56),
            default_km(model, par=km_nadp, rxn=rxn, subs=nadp, value=0.22),
            keq,
        ],
    )

    return model


def add_fnr_mm(
    model: Model,
    *,
    rxn: str | None = None,
    fd_ox: str | None = None,
    fd_red: str | None = None,
    nadph: str | None = None,
    nadp: str | None = None,
    ph: str | None = None,
    kcat: str | None = None,
    km_fd: str | None = None,
    km_nadp: str | None = None,
    e0: str | None = None,
    convf: str | None = None,
) -> Model:
    """Add fnr mm (Ferredoxin-NADP reductase) to model."""
    rxn = default_name(rxn, n.fnr)
    fd_ox = default_name(fd_ox, n.fd_ox)
    fd_red = default_name(fd_red, n.fd_red)
    nadph = default_name(nadph, n.nadph)
    nadp = default_name(nadp, n.nadp)
    ph = default_name(ph, n.ph)

    convf = static(model, n.convf(), 3.2e-2) if convf is None else convf

    model.add_derived(
        n.keq(rxn),
        fn=_keq_fnr,
        args=[
            "E^0_Fd",
            "F",
            "E^0_NADP",
            ph,
            "dG_pH",
            "RT",
        ],
    )

    model.add_reaction(
        name=rxn,
        fn=_rate_fnr_2019,
        stoichiometry={
            fd_ox: 2,
            nadph: Derived(fn=value, args=[convf]),
        },
        args=[
            fd_ox,
            fd_red,
            nadph,
            nadp,
            static(model, n.km(rxn, n.fd_red()), 1.56) if km_fd is None else km_fd,
            static(model, n.km(rxn, n.nadp()), 0.22) if km_nadp is None else km_nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=3.0,  # Source
                kcat_value=500.0,  # Source
            ),
            n.keq(rxn),
            convf,
        ],
    )
    return model


def add_fnr_static(
    model: Model,
    *,
    rxn: str | None = None,
    nadp: str | None = None,
    nadph: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Saadat version to put into Poolman model."""
    rxn = default_name(rxn, n.fnr)
    nadp = default_name(nadp, n.nadp)
    nadph = default_name(nadph, n.nadph)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadp: -1.0,
                nadph: 1.0,
            },
        ),
        args=[
            nadp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=2.816,  # Source
            ),
            default_kms(model, par=kms, rxn=rxn, value=0.19),
        ],
    )

    return model


def add_fnr_energy_dependent(
    model: Model,
    *,
    rxn: str | None = None,
    nadp: str | None = None,
    nadph: str | None = None,
    energy: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add fnr energy dependent (Ferredoxin-NADP reductase) to model."""
    rxn = default_name(rxn, n.fnr)
    nadp = default_name(nadp, n.nadp)
    nadph = default_name(nadph, n.nadph)
    energy = default_name(energy, n.energy)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                # Substrates
                nadp: -1.0,
                energy: -1.0,
                # Products
                nadph: 1.0,
            },
        ),
        args=[
            nadp,
            energy,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=2.816,  # Source
            ),
            default_kms(model, par=kms, rxn=rxn, value=0.19),
        ],
    )

    return model

"""H+-transporting two-sector ATPase.

ADP + Orthophosphate -> ATP

EC 3.6.3.14

Equilibrator
ADP(aq) + Orthophosphate(aq) ⇌ ATP(aq) + H2O(l)
Keq = 6.4e-6 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

import math

from mxlpy import Derived, Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, neg_div, value
from mxlbricks.utils import (
    default_name,
    filter_stoichiometry,
    static,
)


def _keq_atp(
    pH: float,
    DeltaG0_ATP: float,
    dG_pH: float,
    HPR: float,
    pHstroma: float,
    Pi_mol: float,
    RT: float,
) -> float:
    """Equilibrium constant for ATP synthase, driven by the transmembrane proton gradient."""
    delta_g = DeltaG0_ATP - dG_pH * HPR * (pHstroma - pH)
    return Pi_mol * math.exp(-delta_g / RT)


def _rate_atp_synthase_2000(
    adp: float,
    pi: float,
    v16: float,
    km161: float,
    km162: float,
) -> float:
    """ATP synthase rate (Poolman 2000): bi-substrate Michaelis-Menten on ADP and Pi."""
    return v16 * adp * pi / ((adp + km161) * (pi + km162))


def _rate_atp_synthase_2016(
    ATP: float,
    ADP: float,
    Keq_ATPsynthase: float,
    kATPsynth: float,
) -> float:
    """ATP synthase rate (2016 formulation): linear reversible kinetics driven by Keq."""
    return kATPsynth * (ADP - ATP / Keq_ATPsynthase)


def _rate_atp_synthase_2019(
    ATP: float,
    ADP: float,
    Keq_ATPsynthase: float,
    kATPsynth: float,
    convf: float,
) -> float:
    """ATP synthase rate (2019 formulation): same as 2016 but ADP/ATP scaled by convf."""
    return kATPsynth * (ADP / convf - ATP / convf / Keq_ATPsynthase)


def _rate_static_energy(
    adp: float,
    pi: float,
    energy: float,
    v16: float,
    km161: float,
    km162: float,
) -> float:
    """ATP synthase rate scaled by an energy variable representing photosynthetic driving force."""
    return adp * pi * energy * v16 / ((adp + km161) * (pi + km162))


def _atp_gamma(
    Pi: float,
    ATP: float,
    ADP: float,
    convf: float,
) -> float:
    """Mass-action ratio Γ = [ATP] / ([ADP][Pi]) for ATP synthase free-energy calculation."""
    return (ATP / convf) / ((ADP / convf) * (Pi / 1000))


def _delta_g_atp_synthase(
    pH: float,
    gammaATP: float,
    DeltaG0_ATP: float,
    dG_pH: float,
    HPR: float,
    pHstroma: float,
    RT: float,
) -> float:
    """Free-energy change of ATP synthesis including PMF contribution and mass-action correction."""
    return DeltaG0_ATP - dG_pH * HPR * (pHstroma - pH) + RT * math.log(gammaATP)


def _atp_synthase2(
    DeltaGATPsyn: float,
    ATPturnover: float,
) -> float:
    """ATP synthase flux proportional to free-energy drop and turnover number."""
    return -DeltaGATPsyn * ATPturnover


def add_atp_synthase_mmol_chl(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    h_stroma: str | None = None,
    h_lumen: str | None = None,
    ph_stroma: str | None = None,
    ph_lumen: str | None = None,
    kf: str | None = None,
    hpr: str | None = None,
    bh: str | None = None,  # proton buffering
) -> Model:
    """Add atp synthase mmol chl (H+-transporting two-sector ATPase) to model."""
    rxn = default_name(rxn, n.atp_synthase)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    h_stroma = default_name(h_stroma, n.h)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    ph_stroma = default_name(ph_stroma, n.ph)
    ph_lumen = default_name(ph_lumen, lambda: n.ph("_lumen"))

    kf = static(model, n.kf(rxn), 20.0) if kf is None else kf
    hpr = static(model, "HPR", 14.0 / 3.0) if hpr is None else hpr
    bh = static(model, "bH", 100.0) if bh is None else bh

    model.add_parameter("Pi_mol", 0.01)
    model.add_parameter("DeltaG0_ATP", 30.6)

    model.add_derived(
        name=n.keq(rxn),
        fn=_keq_atp,
        args=[
            ph_lumen,
            "DeltaG0_ATP",
            "dG_pH",
            hpr,
            ph_stroma,
            "Pi_mol",
            "RT",
        ],
    )

    model.add_reaction(
        name=rxn,
        fn=_rate_atp_synthase_2016,
        stoichiometry={
            atp: 1.0,
            h_lumen: Derived(
                fn=neg_div,
                args=[
                    "HPR",
                    bh,
                ],
            ),
        },
        args=[
            atp,
            adp,
            n.keq(rxn),
            kf,
        ],
    )
    return model


def add_atp_synthase_mm(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    h_stroma: str | None = None,
    h_lumen: str | None = None,
    ph_stroma: str | None = None,
    ph_lumen: str | None = None,
    kf: str | None = None,
    hpr: str | None = None,
    bh: str | None = None,  # proton buffering
    convf: str | None = None,
) -> Model:
    """Add atp synthase mm (H+-transporting two-sector ATPase) to model."""
    rxn = default_name(rxn, n.atp_synthase)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    h_stroma = default_name(h_stroma, n.h)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    ph_stroma = default_name(ph_stroma, n.ph)
    ph_lumen = default_name(ph_lumen, lambda: n.ph("_lumen"))

    kf = static(model, n.kf(rxn), 20.0) if kf is None else kf
    hpr = static(model, "HPR", 14.0 / 3.0) if hpr is None else hpr
    bh = static(model, "bH", 100.0) if bh is None else bh
    convf = static(model, n.convf(), 3.2e-2) if convf is None else convf

    model.add_parameter("Pi_mol", 0.01)
    model.add_parameter("DeltaG0_ATP", 30.6)

    model.add_derived(
        name=n.keq(rxn),
        fn=_keq_atp,
        args=[
            ph_lumen,
            "DeltaG0_ATP",
            "dG_pH",
            hpr,
            ph_stroma,
            "Pi_mol",
            "RT",
        ],
    )

    model.add_reaction(
        name=rxn,
        fn=_rate_atp_synthase_2019,
        stoichiometry={
            h_lumen: Derived(
                fn=neg_div,
                args=[hpr, bh],
            ),
            atp: Derived(fn=value, args=[convf]),
        },
        args=[
            atp,
            adp,
            n.keq(rxn),
            kf,
            convf,
        ],
    )
    return model


def add_atp_synthase_static_protons(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_adp: str | None = None,
    km_pi: str | None = None,
) -> Model:
    """Add ATP synthase reaction as used in Poolman 2000."""
    rxn = default_name(rxn, n.atp_synthase)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)

    km_adp = (
        static(model, n.km(rxn, adp), 0.014) if km_adp is None else km_adp
    )  # FIXME: source
    km_pi = (
        static(model, n.km(rxn, pi), 0.3) if km_pi is None else km_pi
    )  # FIXME: source
    kcat = static(model, n.kcat(rxn), 2.8) if kcat is None else kcat  # FIXME: source
    e0 = static(model, n.e0(rxn), 1.0) if e0 is None else e0  # FIXME: source
    model.add_derived(vmax := n.vmax(rxn), fn=mass_action_1s, args=[kcat, e0])

    model.add_reaction(
        name=rxn,
        fn=_rate_atp_synthase_2000,
        stoichiometry=filter_stoichiometry(
            model,
            {
                adp: -1.0,
                atp: 1.0,
            },
        ),
        args=[
            adp,
            pi,
            vmax,
            km_adp,
            km_pi,
        ],
    )
    return model


def add_atp_synthase_energy_dependent(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    pi: str | None = None,
    energy: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_adp: str | None = None,
    km_pi: str | None = None,
) -> Model:
    """Add atp synthase energy dependent (H+-transporting two-sector ATPase) to model."""
    rxn = default_name(rxn, n.atp_synthase)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    pi = default_name(pi, n.pi)
    energy = default_name(energy, n.energy)

    km_adp = (
        static(model, n.km(rxn, adp), 0.014) if km_adp is None else km_adp
    )  # FIXME: source
    km_pi = (
        static(model, n.km(rxn, pi), 0.3) if km_pi is None else km_pi
    )  # FIXME: source
    kcat = static(model, n.kcat(rxn), 2.8) if kcat is None else kcat  # FIXME: source
    e0 = static(model, n.e0(rxn), 1.0) if e0 is None else e0  # FIXME: source
    model.add_derived(vmax := n.vmax(rxn), fn=mass_action_1s, args=[kcat, e0])

    model.add_reaction(
        name=rxn,
        fn=_rate_static_energy,
        stoichiometry=filter_stoichiometry(
            model,
            {
                # Substrates
                adp: -1.0,
                energy: -1.0,
                # Products
                atp: 1.0,
            },
        ),
        args=[
            adp,
            pi,
            energy,
            vmax,
            km_adp,
            km_pi,
        ],
    )
    return model


def add_atp_synthase_2024(
    model: Model,
    *,
    rxn: str | None = None,
    atp: str | None = None,
    adp: str | None = None,
    h_lumen: str | None = None,
    ph_stroma: str | None = None,
    ph_lumen: str | None = None,
    pi: str | None = None,
    kf: str | None = None,
    hpr: str | None = None,
    bh: str | None = None,  # proton buffering
    convf: str | None = None,
) -> Model:
    """Add atp synthase 2024 (H+-transporting two-sector ATPase) to model."""
    rxn = default_name(rxn, n.atp_synthase)
    atp = default_name(atp, n.atp)
    adp = default_name(adp, n.adp)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    ph_stroma = default_name(ph_stroma, n.ph)
    ph_lumen = default_name(ph_lumen, lambda: n.ph("_lumen"))
    pi = default_name(pi, n.pi)

    kf = static(model, n.kf(rxn), 20.0) if kf is None else kf
    hpr = static(model, "HPR", 14.0 / 3.0) if hpr is None else hpr
    bh = static(model, "bH", 100.0) if bh is None else bh
    convf = static(model, n.convf(), 3.2e-2) if convf is None else convf

    model.add_parameter(delta_g := "DeltaG0_ATP", 30.6)
    model.add_parameter(atp_turnover := "ATPturnover", 90)

    model.add_derived(
        name="ATP_gamma",
        fn=_atp_gamma,
        args=[
            pi,
            atp,
            adp,
            convf,
        ],
    )

    model.add_derived(
        name="DeltaGATPsyn",
        fn=_delta_g_atp_synthase,
        args=[
            ph_lumen,
            "ATP_gamma",
            delta_g,
            "dG_pH",
            hpr,
            ph_stroma,
            "RT",
        ],
    )

    model.add_reaction(
        name="vATPsynthase",
        fn=_atp_synthase2,
        stoichiometry={
            h_lumen: Derived(fn=neg_div, args=[hpr, bh]),
            atp: Derived(fn=value, args=[convf]),
        },
        args=[
            "DeltaGATPsyn",
            atp_turnover,
        ],
    )
    return model

"""Ribulose-1,5-bisphosphate carboxylase/oxygenase.

Enzyme catalysing both carboxylation as well as oxygenation of ribulose-1,5-bisphosphate
leading to either 2xPGA or 1xPGA and 1xPGO


Equilibrator (carboxylation)
    D-Ribulose 1,5-bisphosphate(aq) + CO2(total) ⇌ 2 3-Phospho-D-glycerate(aq)
    Keq = 1.6e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)

Equilibrator (oxygenation)
    Oxygen(aq) + D-Ribulose 1,5-bisphosphate(aq) ⇌ 3-Phospho-D-glycerate(aq) + 2-Phosphoglycolate(aq)
    Keq = 2.9e91 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)


Following inhibition mechanisms are known
    - PGA (Poolman 2000)
    - FBP (Poolman 2000)
    - SBP (Poolman 2000)
    - Orthophosphate (Poolman 2000)
    - NADPH (Poolman 2000)
    - PGO (FIXME)


Because of it's complex dynamics, multiple kinetic descriptions of rubisco are possible,
some of which have been implemented here.
    - Poolman 2000, doi: FIXME
    - Witzel 2010, doi: FIXME

Kinetic parameters
------------------
kcat (CO2)
    - 3 s^1 (Stitt 2010)

Witzel:
    gamma = 1 / km_co2
    omega = 1 / km_o2
    lr = k_er_minus / k_er_plus
    lc = k_er_minus / (omega * kcat_carb)
    lrc = k_er_minus / (gamma * k_er_plus)
    lro = k_er_minus / (omega * k_er_plus)
    lo = k_er_minus / (omega * k_oxy)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import div, mass_action_1s, mul, one_div
from mxlbricks.utils import (
    default_km,
    default_name,
    default_par,
    default_vmax,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def _rate_poolman_5i(
    rubp: float,
    pga: float,
    co2: float,
    vmax: float,
    kms_rubp: float,
    kms_co2: float,
    # inhibitors
    ki_pga: float,
    fbp: float,
    ki_fbp: float,
    sbp: float,
    ki_sbp: float,
    pi: float,
    ki_p: float,
    nadph: float,
    ki_nadph: float,
) -> float:
    top = vmax * rubp * co2
    btm = (
        rubp
        + kms_rubp
        * (
            1
            + pga / ki_pga
            + fbp / ki_fbp
            + sbp / ki_sbp
            + pi / ki_p
            + nadph / ki_nadph
        )
    ) * (co2 + kms_co2)
    return top / btm


def _rate_witzel_5i(
    rubp: float,
    s2: float,
    vmax: float,
    gamma_or_omega: float,
    co2: float,
    o2: float,
    lr: float,
    lc: float,
    lo: float,
    lrc: float,
    lro: float,
    i1: float,  # pga
    ki1: float,
    i2: float,  # fbp
    ki2: float,
    i3: float,  # sbp
    ki3: float,
    i4: float,  # pi
    ki4: float,
    i5: float,  # nadph
    ki5: float,
) -> float:
    vmax_app = (gamma_or_omega * vmax * s2 / lr) / (1 / lr + co2 / lrc + o2 / lro)
    km_app = 1 / (1 / lr + co2 / lrc + o2 / lro)
    return (vmax_app * rubp) / (
        rubp
        + km_app
        * (
            1
            + co2 / lc
            + o2 / lo
            + i1 / ki1
            + i2 / ki2
            + i3 / ki3
            + i4 / ki4
            + i5 / ki5
        )
    )


def add_rubisco_poolman(
    model: Model,
    *,
    rxn: str | None = None,
    rubp: str | None = None,
    pga: str | None = None,
    co2: str | None = None,
    fbp: str | None = None,
    sbp: str | None = None,
    pi: str | None = None,
    nadph: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_co2: str | None = None,
    km_rubp: str | None = None,
    ki_pga: str | None = None,
    ki_fbp: str | None = None,
    ki_sbp: str | None = None,
    ki_pi: str | None = None,
    ki_nadph: str | None = None,
) -> Model:
    """Add rubisco poolman (Ribulose-1,5-bisphosphate carboxylase/oxygenase) to model."""
    rxn = default_name(rxn, n.rubisco_carboxylase)
    rubp = default_name(rubp, n.rubp)
    pga = default_name(pga, n.pga)
    co2 = default_name(co2, n.co2)
    fbp = default_name(fbp, n.fbp)
    sbp = default_name(sbp, n.sbp)
    pi = default_name(pi, n.pi)
    nadph = default_name(nadph, n.nadph)

    model.add_reaction(
        name=n.rubisco_carboxylase(),
        fn=_rate_poolman_5i,
        stoichiometry=filter_stoichiometry(
            model,
            {
                rubp: -1.0,
                pga: 2.0,
                co2: -1,
            },
        ),
        args=[
            rubp,
            pga,
            co2,
            default_vmax(
                model, e0=e0, kcat=kcat, rxn=rxn, e0_value=1.0, kcat_value=0.34 * 8
            ),
            default_km(model, par=km_rubp, rxn=rxn, subs=rubp, value=0.02),
            default_km(model, par=km_co2, rxn=rxn, subs=co2, value=0.0107),
            default_par(model, par=ki_pga, name=n.ki(rxn, n.pga()), value=0.04),
            fbp,
            default_par(model, par=ki_fbp, name=n.ki(rxn, n.fbp()), value=0.04),
            sbp,
            default_par(model, par=ki_sbp, name=n.ki(rxn, n.sbp()), value=0.075),
            pi,
            default_par(model, par=ki_pi, name=n.ki(rxn, n.pi()), value=0.9),
            nadph,
            default_par(model, par=ki_nadph, name=n.ki(rxn, n.nadph()), value=0.07),
        ],
    )
    return model


def add_rubisco(
    model: Model,
    *,
    rxn_carb: str | None = None,
    rxn_ox: str | None = None,
    rubp: str | None = None,
    co2: str | None = None,
    o2: str | None = None,
    pga: str | None = None,
    pgo: str | None = None,
    fbp: str | None = None,
    sbp: str | None = None,
    pi: str | None = None,
    nadph: str | None = None,
    kcat_carb: str | None = None,
    kcat_ox: str | None = None,
    e0: str | None = None,
    km_co2: str | None = None,
    km_o2: str | None = None,
    km_rubp: str | None = None,
    k_er_plus: str | None = None,
    k_er_minus: str | None = None,
    ki_pga: str | None = None,
    ki_fbp: str | None = None,
    ki_sbp: str | None = None,
    ki_pi: str | None = None,
    ki_nadph: str | None = None,
) -> Model:
    """Add rubisco (Ribulose-1,5-bisphosphate carboxylase/oxygenase) to model."""
    enzyme = n.rubisco()
    rxn_carb = default_name(rxn_carb, n.rubisco_carboxylase)
    rxn_ox = default_name(rxn_ox, n.rubisco_oxygenase)
    rubp = default_name(rubp, n.rubp)
    co2 = default_name(co2, n.co2)
    o2 = default_name(o2, n.o2)
    pga = default_name(pga, n.pga)
    pgo = default_name(pgo, n.pgo)
    fbp = default_name(fbp, n.fbp)
    sbp = default_name(sbp, n.sbp)
    pi = default_name(pi, n.pi)
    nadph = default_name(nadph, n.nadph)

    kcat_carb = (
        static(model, n.kcat(rxn_carb), 3.1) if kcat_carb is None else kcat_carb
    )  # FIXME: source
    kcat_ox = (
        static(model, n.kcat(rxn_ox), 1.125) if kcat_carb is None else kcat_carb
    )  # FIXME: source
    km_co2 = (
        static(model, n.km(enzyme, co2), 10.7 / 1000) if km_co2 is None else km_co2
    )  # FIXME: source
    km_o2 = (
        static(model, n.km(enzyme, o2), 295 / 1000) if km_o2 is None else km_o2
    )  # FIXME: source
    km_rubp = (
        static(model, n.km(enzyme, rubp), 0.02) if km_rubp is None else km_rubp
    )  # FIXME: source
    e0 = static(model, n.e0(enzyme), 0.16) if e0 is None else e0  # FIXME: source

    model.add_derived(
        vmax_carb := n.vmax(rxn_carb), fn=mass_action_1s, args=[kcat_carb, e0]
    )
    model.add_derived(vmax_ox := n.vmax(rxn_ox), fn=mass_action_1s, args=[kcat_ox, e0])

    ki_pga = static(model, n.ki(enzyme, n.pga()), 0.04) if ki_pga is None else ki_pga
    ki_fbp = static(model, n.ki(enzyme, n.fbp()), 0.04) if ki_fbp is None else ki_fbp
    ki_sbp = static(model, n.ki(enzyme, n.sbp()), 0.075) if ki_sbp is None else ki_sbp
    ki_pi = static(model, n.ki(enzyme, n.pi()), 0.9) if ki_pi is None else ki_pi
    ki_nadph = (
        static(model, n.ki(enzyme, n.nadph()), 0.07) if ki_nadph is None else ki_nadph
    )

    k_er_plus = (
        static(model, "k_er_plus", 0.15 * 1000) if k_er_plus is None else k_er_plus
    )  # 1 / (mM * s)
    k_er_minus = (
        static(model, "k_er_minus", 0.0048) if k_er_minus is None else k_er_minus
    )  # 1 / s

    model.add_derived(gamma := "gamma", one_div, args=[km_co2])
    model.add_derived(omega := "omega", one_div, args=[km_o2])
    model.add_derived(
        omega_kcat_carb := "omega_kcat_carb",
        mul,
        args=[omega, n.kcat(n.rubisco_carboxylase())],
    )
    model.add_derived(
        omega_koxy := "omega_koxy", mul, args=[omega, n.kcat(n.rubisco_oxygenase())]
    )
    model.add_derived(omega_ker_plus := "omega_ker_plus", mul, args=[omega, k_er_plus])
    model.add_derived(gamma_ker_plus := "gamma_ker_plus", mul, args=[gamma, k_er_plus])
    model.add_derived(lr := "lr", div, args=[k_er_minus, k_er_plus])
    model.add_derived(lc := "lc", div, args=[k_er_minus, omega_kcat_carb])
    model.add_derived(lrc := "lrc", div, args=[k_er_minus, gamma_ker_plus])
    model.add_derived(lro := "lro", div, args=[k_er_minus, omega_ker_plus])
    model.add_derived(lo := "lo", div, args=[k_er_minus, omega_koxy])
    model.add_reaction(
        name=n.rubisco_carboxylase(),
        fn=_rate_witzel_5i,
        stoichiometry=filter_stoichiometry(
            model,
            {
                rubp: -1.0,
                pga: 2.0,
                co2: -1,
            },
        ),
        args=[
            rubp,
            co2,
            vmax_carb,
            gamma,  # 1 / km_co2
            co2,
            o2,
            lr,
            lc,
            lo,
            lrc,
            lro,
            pga,
            ki_pga,
            fbp,
            ki_fbp,
            sbp,
            ki_sbp,
            pi,
            ki_pi,
            nadph,
            ki_nadph,
        ],
    )
    model.add_reaction(
        name=n.rubisco_oxygenase(),
        fn=_rate_witzel_5i,
        stoichiometry=filter_stoichiometry(
            model,
            {
                rubp: -1.0,
                o2: -1.0,
                pga: 1.0,
                pgo: 1.0,
            },
        ),
        args=[
            rubp,
            o2,
            vmax_ox,
            omega,  # 1 / km_o2
            co2,
            o2,
            lr,
            lc,
            lo,
            lrc,
            lro,
            pga,
            ki_pga,
            fbp,
            ki_fbp,
            sbp,
            ki_sbp,
            pi,
            ki_pi,
            nadph,
            ki_nadph,
        ],
    )

    return model

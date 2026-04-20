from typing import cast

import numpy as np
from mxlpy import Derived, InitialAssignment, Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, mass_action_2s, moiety_1, value
from mxlbricks.models import get_saadat2021


### PSII and PSI ###
def _moiety_3(
    c1: float,
    c2: float,
    c3: float,
    total: float,
) -> float:
    return total - c1 - c2 - c3


def _normalize_concentration(
    concentration: float,
    total: float,
) -> float:
    return concentration / total


def _normalize_2_concentrations(
    c1: float,
    c2: float,
    total: float,
) -> float:
    return (c1 + c2) / total


def _mass_action_22_rev(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    kf: float,
    keq: float,
) -> float:
    return kf * s1 * s2 - (kf / keq) * p1 * p2


def _kquencher(
    s: float,
    q: float,
    kH_Qslope: float,
    kH0: float,
) -> float:
    return (kH0 + kH_Qslope * q) * s


def _v_ps1(
    P700FA: float,
    ps2cs: float,
    pfd: float,
) -> float:
    return (1 - ps2cs) * pfd * P700FA


def _v_mehler(
    PSI_red_acceptor: float,
    O2ext: float,
    kMehler: float,
) -> float:
    return kMehler * O2ext * PSI_red_acceptor


def _fluo(
    Q: float,
    B0: float,
    B2: float,
    ps2cs: float,
    k2: float,
    kF: float,
    kH_Qslope: float,
    kH0: float,
) -> float:
    kH = kH0 + kH_Qslope * Q
    return (ps2cs * kF * B0) / (kF + k2 + kH) + (ps2cs * kF * B2) / (kF + kH)


def _keq_pcp700(
    E0_PC: float,
    F: float,
    E0_P700: float,
    RT: float,
) -> float:
    DG = -(-E0_PC * F) + (-E0_P700 * F)
    return np.exp(-DG / RT)


def _keq_faf_d(
    E0_FA: float,
    F: float,
    E0_Fd: float,
    RT: float,
) -> float:
    DG = -(-E0_FA * F) + (-E0_Fd * F)
    return np.exp(-DG / RT)


### b6f ###


def _four_div_by(
    x: float,
) -> float:
    return 4.0 / x


def _k_b6f(
    pH: float,
    pKreg: float,
    b6f_content: float,
    max_b6f: float,
) -> float:
    pHmod = 1 - (1 / (10 ** (pH - pKreg) + 1))
    return pHmod * b6f_content * max_b6f


def _vb6f_2024(
    PC: float,
    PCred: float,
    PQ: float,
    PQred: float,
    _k_b6f: float,
    Keq: float,
) -> float:
    f = PQred / (PQred + PQ)
    return f * PC * _k_b6f - (1 - f) * PCred * (_k_b6f / Keq)


def _reversible_mass_action_1s_1p(
    s: float,
    p: float,
    kf: float,
    kb: float,
) -> float:
    return kf * s - kb * p


def _half(
    x: float,
) -> float:
    return x / 2


def _protons_lumen(
    pH_lumen: float,
) -> float:
    return (10 ** (-pH_lumen)) / 2.5e-4


def _protons_stroma(
    pH_stroma: float,
) -> float:
    return (10 ** (-pH_stroma)) / 3.2e-5


def _v_at_psynthase_mod(
    ATP: float,
    ATP_activity: float,
    _atp_pmf_activity: float,
    k_ATPsynthase: float,
    ADP: float,
    _keq_atp: float,
    convf: float,
) -> float:
    return (
        ATP_activity
        * _atp_pmf_activity
        * k_ATPsynthase
        * (ADP / convf - ATP / convf / _keq_atp)
    )


def _atp_pmf_activity(
    pK0E: float,
    b: float,
    pH_lumen: float,
    pH: float,
    F: float,
    RT: float,
    delta_psi: float,
) -> float:
    _pmf = delta_psi - np.log(10) * ((RT) / F) * (pH_lumen - pH)
    x = np.log(10 ** (-pK0E)) + b * (_pmf * F) / (RT)
    return (np.e**x) / (1 + np.e**x)


def _v_at_pactivity(
    ATPactivity: float,
    light: float,
    kActATPase: float,
    kDeactATPase: float,
) -> float:
    """Activation of ATPsynthase by light"""
    if light > 0.0:
        return kActATPase * (1 - ATPactivity)
    return -kDeactATPase * ATPactivity


def _atp_pmf_activity2(
    pK0E: float,
    b: float,
    pH_lumen: float,
    pH: float,
    F: float,
    RT: float,
    delta_psi: float,
) -> float:
    _pmf = delta_psi - np.log(10) * ((RT) / F) * (pH_lumen - pH)
    x = np.log(10 ** (-pK0E)) + b * (_pmf * F) / (RT)
    return (np.e**x) / (1 + np.e**x)


def _deltap_h(
    pH: float,
    pH_lumen: float,
    dG: float,
) -> float:
    return dG * (pH - pH_lumen)


def _pmf(
    _deltap_h: float,
    delta_psi: float,
    F: float,
) -> float:
    return F * delta_psi + _deltap_h


def _pmf_in_v(
    delta_psi: float,
    pH_lumen: float,
    pH: float,
    RT: float,
    F: float,
) -> float:
    return delta_psi - np.log(10) * ((RT) / F) * (pH_lumen - pH)


def _voltage_turnover_mol_chl_per_mmol(
    capacitance_specific: float,
    molChl_per_area_membrane: float,
    F: float,
) -> float:
    area_permolChl = 1 / molChl_per_area_membrane
    return F / (capacitance_specific * area_permolChl)


def _initial_delta_psi(
    pH: float,
    pH_lumen: float,
    R: float,
    F: float,
    T: float,
) -> float:
    """
    Estimation of delta psi in the dark - assuming delta_pH and delta_psi have equal contribution to pmf
    """
    return np.log(10) * ((R * T) / F) * (pH - pH_lumen)


def _keq_atp(
    _pmf: float,
    DeltaG0_ATP: float,
    Pi_mol: float,
    HPR: float,
    RT: float,
) -> float:
    DG = DeltaG0_ATP - HPR * _pmf

    return Pi_mol * np.exp(-DG / RT)


def _keq_cytb6f(
    pH: float,
    _pmf: float,
    F: float,
    E0_PQ: float,
    E0_PC: float,
    RT: float,
    dG_pH: float,
) -> float:
    DG1 = -2 * F * E0_PQ
    DG2 = -F * E0_PC
    DG = -(DG1 + 2 * dG_pH * pH) + 2 * DG2 + 2 * _pmf
    return np.exp(-DG / RT)


def _one_div(
    x: float,
) -> float:
    return x


def _two_div(
    x: float,
) -> float:
    return 2 * x


def _three_div(
    x: float,
) -> float:
    return 3 * x


def _neg_one_div(
    x: float,
) -> float:
    return -1 * x


def _atp_div(
    HPR: float,
    x: float,
) -> float:
    return -HPR * x


def _four_div(
    x: float,
) -> float:
    return 4 * x


def _reg_kea(
    pH: float,
    ATP: float,
    KEA3_pH_reg: float,
    KEA3_ATP_treshold: float,
) -> float:
    pH_inhib = (1 - 0.1) / (1 + np.exp((pH - KEA3_pH_reg) / 0.001))
    ATP_inhib = (1 - 0.1) / (1 + np.exp((KEA3_ATP_treshold - ATP) / 0.01))
    return pH_inhib * ATP_inhib


def _dg_k(
    Klumen: float,
    Kstroma: float,
    delta_psi: float,
    RT: float,
    F: float,
) -> float:
    return (-(RT / F) * np.log10(Kstroma / Klumen) + delta_psi) * F


def _cl_driving_force(
    delta_psi: float,
    Cl_lumen: float,
    Cl_stroma: float,
    RT: float,
    F: float,
) -> float:
    return ((RT / F) * np.log10(Cl_stroma / Cl_lumen) + delta_psi) * F


def _keq_ndh1(
    _pmf: float,
    E0_Fd: float,
    F: float,
    E0_PQ: float,
    pHstroma: float,
    dG_pH: float,
    RT: float,
) -> float:
    DG1 = -E0_Fd * F
    DG2 = -2 * E0_PQ * F
    DG = -2 * DG1 + DG2 + 2 * dG_pH * pHstroma + 4 * _pmf
    return np.exp(-DG / RT)


def _v_kea(
    Klumen: float,
    H: float,
    Kstroma: float,
    k_KEA: float,
    Hstroma: float,
    _reg_kea: float,
) -> float:
    v_KEA = k_KEA * (H * Kstroma - Hstroma * Klumen) * _reg_kea
    return max(
        v_KEA,
        0,
    )


def _v_voltage_k_channel(
    delta_psi_ions: float,
    Klumen: float,
    Kstroma: float,
    _dg_k: float,
    perm_K: float,
    K_delta_psi_treshold: float,
) -> float:
    voltage_dependence = (1 - 0.1) / (
        1 + np.exp(-(delta_psi_ions - K_delta_psi_treshold) / 0.001)
    )
    return (
        perm_K * _dg_k * voltage_dependence * (Klumen / Kstroma)
    )  # why divided K_total/2


def _v_vccn1(
    Cl_stroma: float,
    Cl_lumen: float,
    _cl_driving_force: float,
    delta_psi_ions: float,
    k_VCCN1: float,
    VCCN_delta_psi_treshold: float,
) -> float:
    voltage_gate = (1 - 0.1) / (
        1 + np.exp(-(delta_psi_ions - VCCN_delta_psi_treshold) / 0.001)
    )
    return voltage_gate * k_VCCN1 * _cl_driving_force * (Cl_stroma / Cl_lumen)


def _v_cl_ce(
    Cl_lumen: float,
    Cl_stroma: float,
    kClCe: float,
    PQ: float,
    _cl_driving_force: float,
    ClCe_PQ: float,
) -> float:  # correct
    activation = (0.2 - 0.1) / (1 + np.exp(-(PQ - ClCe_PQ) / 0.1))
    return kClCe * activation * _cl_driving_force * (Cl_stroma / Cl_lumen)


def _cl_ce_activation(
    ATP: float,
    ClCe_ATP_threshold: float,
) -> float:
    return (1 - 0.1) / (1 + np.exp((ATP - ClCe_ATP_threshold) / 0.01))


def _cl_ce_bi(
    Cl_lumen: float,
    Cl_stroma: float,
    kClCe: float,
    activation: float,
) -> float:  # correct
    return kClCe * (Cl_stroma - Cl_lumen) * activation


def _cl_ce_ch(
    Cl_lumen: float,
    Cl_stroma: float,
    kClCe: float,
    activation: float,
    protons: float,
    _protons_lumen: float,
) -> float:  # correct
    return kClCe * ((Cl_lumen * protons) / (_protons_lumen * Cl_stroma)) * activation


def _v_cl_leak(
    kCl_leak: float,
    Cl_lumen: float,
    Cl_stroma: float,
    PQ: float,
    Cl_leak_PQ: float,
    total_div: float,
) -> float:
    activation = (1 - 0.1) / (1 + np.exp(-(PQ - Cl_leak_PQ) / 0.1))
    return kCl_leak * ((Cl_lumen - Cl_stroma) ** 2) / (total_div) * activation


def _v_ndh1(
    A1: float,
    Fdred: float,
    PQ: float,
    pHlumen: float,
    kNDH1: float,
) -> float:
    return (
        kNDH1
        * ((Fdred**2) * PQ)
        * ((1 - 0.1) / (1 + np.exp(-((A1 - 0.02) / 0.01))))
        * (10 ** (pHlumen - 6.5) / (10 ** (pHlumen - 6.5) + 0.5))
    )


def _squared(
    x: float,
) -> float:
    return x**2


def _build_full_model(
    m: Model,
    *,
    chl_lumen: str = "_lumen",
    lumen_reactions: list[str] | None = None,
    stroma_reactions: list[str] | None = None,
    ClCe: str = "exporter",
) -> Model:
    """
    Changes:
    1. PSII + PSI dynamics
    2. b6f module photosynthetic control
    3. ATP synthase
    4. tracking pH and variable stroma
    3. delta psi and ion channels
    """

    lumen_reactions = (
        ["B12", "b6f", "proton_leak", "atp_synthase"]
        if lumen_reactions is None
        else lumen_reactions
    )
    stroma_reactions = (
        ["B20", "b6f", "atp_synthase", "proton_leak"]
        if stroma_reactions is None
        else stroma_reactions
    )

    ### PSII + PSI dynamics ###
    m.remove_derived("keq_PCP700")
    if "ps2states" in m._ids:  # noqa: SLF001
        m.remove_surrogate("ps2states")
    m.remove_surrogate("ps1states")
    m.remove_reaction(n.ps2())
    m.remove_reaction(n.ps1())
    m.remove_reaction(n.mehler())
    m.remove_reaction(n.ferredoxin_reductase())
    m.remove_derived("vmax_ferredoxin_reductase")
    m.remove_derived("keq_ferredoxin_reductase")
    m.remove_readout("Fluo")
    m.remove_parameter("kcat_b6f")
    m.remove_parameter("E0_ferredoxin_reductase")
    m.remove_parameter("kcat_ferredoxin_reductase")
    # m.add_parameter("PSII_total", 2.5)
    # m.add_parameter("PSI_total", 2.5)
    m.add_parameter(
        "kH_Qslope",
        5e9,
    )
    m.add_variables(
        {
            "P700FA": 1.506615384275408,  # eq at pfd 800       #"PSItot": 2.5, (in parameter vector of Matuszynska)
            "P700+FA-": 0.019197449388051676,
            "P700FA-": 0.028144516332212766,
            n.b0(): 1.9379789566530539,  # eq at pfd 800
            n.b1(): 9.786232812526368e-08,
            n.b2(): 0.5620208537555176,
        }
    )
    m.add_derived(
        n.keq("PCP700"),
        _keq_pcp700,
        args=["E^0_PC", "F", "E^0_P700", "RT"],
    )
    m.add_derived(
        n.keq("FAFd"),
        _keq_faf_d,
        args=["E^0_FA", "F", "E^0_Fd", "RT"],
    )
    m.add_derived(
        n.b3(),
        _moiety_3,
        args=[n.b0(), n.b1(), n.b2(), "PSII_total"],
    )
    m.add_derived(
        "P700+FA",
        _moiety_3,
        args=["P700FA-", "P700FA", "P700+FA-", "PSI_total"],
    )
    m.add_derived(
        "rel_P700+FA",
        _normalize_concentration,
        args=["P700+FA", "PSI_total"],
    )
    m.add_derived(
        "rel_P700FA",
        _normalize_concentration,
        args=["P700FA", "PSI_total"],
    )
    m.add_derived(
        "rel_P700FA-",
        _normalize_concentration,
        args=["P700FA-", "PSI_total"],
    )
    m.add_derived(
        "rel_P700+FA-",
        _normalize_concentration,
        args=["P700+FA-", "PSI_total"],
    )
    m.add_derived(
        "rel_P700",
        _normalize_2_concentrations,
        args=["P700+FA-", "P700+FA", "PSI_total"],
    )
    m.add_derived(
        "rel_P700+",
        _normalize_2_concentrations,
        args=["P700+FA-", "P700+FA", "PSI_total"],
    )
    m.add_derived(
        "rel_B0",
        _normalize_concentration,
        args=[n.b0(), "PSII_total"],
    )
    m.add_derived(
        "rel_B1",
        _normalize_concentration,
        args=[n.b1(), "PSII_total"],
    )
    m.add_derived(
        "rel_B2",
        _normalize_concentration,
        args=[n.b2(), "PSII_total"],
    )
    m.add_derived(
        "rel_B3",
        _normalize_concentration,
        args=[n.b3(), "PSII_total"],
    )
    m.add_derived(
        n.fluorescence(),
        _fluo,
        args=[n.quencher(), n.b0(), n.b2(), n.ps2cs(), "k2", "kF", "kH_Qslope", "kH0"],
    )

    m.add_reaction(
        "toP700FA-",
        _mass_action_22_rev,
        stoichiometry={"P700+FA-": -1, "P700FA-": 1, n.pc_ox(): 1},
        args=["P700+FA-", n.pc_red(), n.pc_ox(), "P700FA-", "kPCox", n.keq("PCP700")],
    )

    m.add_reaction(
        "toP700FA_v3",
        _mass_action_22_rev,
        stoichiometry={"P700FA-": -1, n.fd_ox(): -1, "P700FA": 1},
        args=["P700FA-", n.fd_ox(), "P700FA", n.fd_red(), "kFdred", n.keq("FAFd")],
    )

    m.add_reaction(
        "toP700+FA",
        _mass_action_22_rev,
        stoichiometry={"P700+FA-": -1, n.fd_ox(): -1},
        args=["P700+FA-", n.fd_ox(), "P700+FA", n.fd_red(), "kFdred", n.keq("FAFd")],
    )

    m.add_reaction(
        "toP700FA_v5",
        _mass_action_22_rev,
        stoichiometry={"P700FA": 1, n.pc_ox(): 1},
        args=["P700+FA", n.pc_red(), "P700FA", n.pc_ox(), "kPCox", n.keq("PCP700")],
    )

    m.add_reaction(
        "PSI",
        _v_ps1,
        stoichiometry={"P700FA": -1, "P700+FA-": 1},
        args=["P700FA", n.ps2cs(), n.pfd()],
    )

    m.add_reaction(
        "mehler1",
        _v_mehler,
        stoichiometry={
            n.h2o2(): Derived(
                fn=value,
                args=["convf"],
            ),
            "P700FA": 2,
            "P700FA-": -2,
        },
        args=["P700FA-", n.o2(chl_lumen), "kMehler"],
    )

    m.add_reaction(
        "mehler2",
        _v_mehler,
        stoichiometry={
            n.h2o2(): Derived(
                fn=value,
                args=["convf"],
            ),
            "P700+FA-": -2,
        },
        args=["P700+FA-", n.o2(chl_lumen), "kMehler"],
    )

    m.add_reaction(
        "B01",
        mass_action_2s,
        stoichiometry={n.b0(): -1, n.b1(): 1},
        args=[n.b0(), n.ps2cs(), n.pfd()],
    )

    m.add_reaction(
        "B10Q",
        _kquencher,
        stoichiometry={n.b1(): -1, n.b0(): 1},
        args=[n.b1(), n.quencher(), "kH_Qslope", "kH0"],
    )

    m.add_reaction(
        "B10F",
        mass_action_1s,
        stoichiometry={n.b1(): -1, n.b0(): 1},
        args=[n.b1(), "kF"],
    )

    m.add_reaction(
        "B12",
        mass_action_1s,
        stoichiometry={
            n.b1(): -1,
            n.b2(): 1,
            n.h(chl_lumen): Derived(
                fn=_one_div,
                args=["bH"],
            ),
        },  # watersplitting occurs on lumenal side and protons will be buffered in the lumen by amino acids
        args=[n.b1(), "k2"],
    )

    m.add_reaction(
        "B20",
        _mass_action_22_rev,
        stoichiometry={n.b2(): -1, n.pq_ox(): -0.5, n.b0(): 1},  # PQred is alm
        args=[n.b2(), n.pq_ox(), n.pq_red(), n.b0(), "kPQred", n.keq(n.pq_red())],
    )

    m.add_reaction(
        "B23",
        mass_action_2s,
        stoichiometry={n.b2(): -1},  # "B3" is alm
        args=[n.b2(), n.ps2cs(), n.pfd()],  # kLII = ps2cs * pfd  and  vB23 = B2 * kLII
    )
    m.add_reaction(
        "B32F",
        mass_action_1s,
        stoichiometry={n.b2(): 1},  # "B3": -1 (alm)
        args=[n.b3(), "kF"],
    )

    m.add_reaction(
        "B32Q",
        _kquencher,
        stoichiometry={n.b2(): 1},  # "B3": -1 (alm)
        args=[n.b3(), n.quencher(), "kH_Qslope", "kH0"],
    )

    ### b6f ###
    bH = m.get_raw_parameters()["bH"].value

    m.add_parameter(
        "b6f_content",
        1,
    )
    m.add_parameter(
        "max_b6f",
        500,
    )
    m.add_parameter(
        "pKreg",
        6.5,
    )
    m.remove_reaction("b6f")

    """m.add_derived("keq_b6f",fn=_keq_cytb6f,
                    args=["pH_lumen", "F", "E^0_PQ", "E^0_PC", "pH", "RT", "dG_pH"])"""

    m.add_derived(
        "keq_b6f_dyn",
        fn=_k_b6f,
        args=[n.ph(chl_lumen), "pKreg", "b6f_content", "max_b6f"],
    )

    m.add_reaction(
        "b6f",
        fn=_vb6f_2024,
        stoichiometry={
            n.pc_ox(): -2,
            n.pq_ox(): 1,
            n.h(chl_lumen): Derived(
                fn=_four_div_by,
                args=["bh"],
            ),
        },
        args=[n.pc_ox(), n.pc_red(), n.pq_ox(), n.pq_red(), "keq_b6f_dyn", "keq_b6f"],
    )

    ### tracking lumen pH ###

    m.remove_derived("pH_lumen")
    m.remove_variable("protons_lumen")
    m.add_variable(
        "pH_lumen",
        6.5,
    )
    m.add_derived(
        "protons_lumen",
        _protons_lumen,
        args=["pH_lumen"],
    )

    HPR = cast(
        float,
        m.get_raw_parameters()["HPR"].value,
    )

    for r in lumen_reactions:
        stoich = m.get_raw_reactions()[r].stoichiometry  # type: ignore

        if r == "B12":
            stoich["pH_lumen"] = -1 / bH  # type: ignore
        elif r == "b6f":
            stoich["pH_lumen"] = -4 / bH  # type: ignore
        elif r == "proton_leak":
            stoich["pH_lumen"] = 1 / bH  # type: ignore
        elif r == "atp_synthase":
            stoich["pH_lumen"] = HPR / bH  # type: ignore
        else:
            continue
        m.update_reaction(
            r,
            stoichiometry=stoich,
        )

    ### tracking stroma pH + variability ###
    m.remove_parameter("protons")
    m.remove_parameter("pH")
    m.add_parameter(
        "stroma_buffering",
        400,
    )
    m.add_variable(
        "pH",
        7.8,
    )
    m.add_derived(
        "protons",
        _protons_stroma,
        args=["pH"],
    )

    buffering = cast(
        float,
        m.get_raw_parameters()["stroma_buffering"].value,
    )

    for r in stroma_reactions:
        stoich = cast(
            dict,
            m.get_raw_reactions()[r].stoichiometry,
        )

        if r == "B20":
            stoich["pH"] = 1 / buffering
        elif r == "b6f":
            stoich["pH"] = 4 / buffering
        elif r == "atp_synthase":
            stoich["pH"] = -HPR / buffering
        elif r == "proton_leak":
            stoich["pH"] = -1 / buffering
        else:
            continue

        m.update_reaction(
            r,
            stoichiometry=stoich,
        )
    ### ATP changes ###

    m.add_variable(
        "ATPactivity",
        0,
    )
    m.add_parameters(
        {
            "kActATPase": 0.01,
            "kDeactATPase": 0.002,
            "k_ATPsynthase": 20,
            "b": 1.7,
            "pK0E": 5.9,
        }
    )

    m.add_reaction(
        "vATPactivity",
        _v_at_pactivity,
        args=["ATPactivity", n.pfd(), "kActATPase", "kDeactATPase"],
        stoichiometry={"ATPactivity": 1},
    )

    stoich = m.get_raw_reactions()["atp_synthase"].stoichiometry  # type: ignore
    m.update_reaction(
        "atp_synthase",
        _v_at_psynthase_mod,
        stoichiometry=stoich,
        args=[
            "ATP",
            "ATPactivity",
            "ATP_pmf_activity",
            "k_ATPsynthase",
            "ADP",
            "kf_atp_synthase",
            "convf",
        ],
    )

    m.add_derived(
        "ATP_pmf_activity",
        _atp_pmf_activity,
        args=["pK0E", "b", "pH_lumen", "pH", "F", "RT", "delta_psi"],
    )

    m.remove_reaction("ex_atp")
    m.remove_reaction("ex_nadph")

    m.update_parameters(
        {
            n.kf(n.ex_atp()): 0.5,
            n.kf(n.ex_nadph()): 0.5,
        }
    )

    m.add_parameters({"k_import_ATP": 0.5, "k_import_NADPH": 0.5})

    m.add_reaction(
        "vATP_shuttle",
        _reversible_mass_action_1s_1p,
        stoichiometry={"ATP": 1},
        args=["ADP", "ATP", "k_import_ATP", n.kf(n.ex_atp())],
    )

    m.add_reaction(
        "vNADPH_shuttle",
        _reversible_mass_action_1s_1p,
        stoichiometry={"NADPH": 1},
        args=["NADP", "NADPH", "k_import_NADPH", n.kf(n.ex_nadph())],
    )

    ### delta psi and ion channels ###
    m.add_parameters(
        {
            "volts_per_charge": 0.000769481926574965,
            "ClCe_PQ": 15,
            "Cl_leak_PQ": 15,
            "KEA3_ATP_treshold": 0.4,
            "KEA3_pH_reg": 7.86,
            "K_delta_psi_treshold": 0.1,
            "VCCN_delta_psi_treshold": 0.08,
            "k_Cl_leak": 0.1,
            "k_NDH1": 1,
            "k_KEA": 200,
            "perm_K": 1,
            "k_VCCN1": 2,
            "k_ClCe": 0.1,
            "K_total": 100,
            "Cl_total": 80,
        }
    )

    m.add_derived(
        "deltapH",
        _deltap_h,
        args=["pH", "pH_lumen", "dG_pH"],
    )

    m.add_derived(
        "deltapH_in_volts",
        _initial_delta_psi,
        args=["pH", "pH_lumen", "R", "F", "T"],
    )

    m.add_variable(
        "delta_psi",
        initial_value=InitialAssignment(
            fn=_initial_delta_psi,
            args=["pH", "pH_lumen", "R", "F", "T"],
        ),
    )

    m.add_derived(
        "pmf",
        _pmf,
        args=["deltapH", "delta_psi", "F"],
    )

    m.add_derived(
        "pmf_in_V",
        _pmf_in_v,
        args=["delta_psi", "pH_lumen", "pH", "RT", "F"],
    )

    m.update_derived(
        "ATP_pmf_activity",
        _atp_pmf_activity2,
        args=["pK0E", "b", "pH_lumen", "pH", "F", "RT", "delta_psi"],
    )

    m.remove_derived("keq_b6f")

    m.add_derived(
        "keq_b6f",
        _keq_cytb6f,
        args=["pH_lumen", "pmf_in_V", "F", "E^0_PQ", "E^0_PC", "RT", "dG_pH"],
    )

    updates = {
        "B12": {
            "delta_psi": Derived(
                fn=_one_div,
                args=["volts_per_charge"],
            )
        },
        "b6f": {
            "delta_psi": Derived(
                fn=_four_div,
                args=["volts_per_charge"],
            )
        },
        "atp_synthase": {
            "delta_psi": Derived(
                fn=_atp_div,
                args=["HPR", "volts_per_charge"],
            )
        },
        "proton_leak": {
            "delta_psi": Derived(
                fn=_neg_one_div,
                args=["volts_per_charge"],
            )
        },
    }

    for r, delta_update in updates.items():
        stoich: dict[str, float | Derived] = cast(
            dict,
            m.get_raw_reactions()[r].stoichiometry,
        )
        stoich.update(delta_update)
        m.update_reaction(
            r,
            stoichiometry=stoich,
        )

    m.add_variable(
        "K_stroma",
        initial_value=InitialAssignment(
            fn=_half,
            args=["K_total"],
        ),
    )
    m.add_variable(
        "Cl_stroma",
        initial_value=InitialAssignment(
            fn=_half,
            args=["Cl_total"],
        ),
    )

    m.add_derived(
        "K_lumen",
        moiety_1,
        args=["K_stroma", "K_total"],
    )
    m.add_derived(
        "Cl_lumen",
        moiety_1,
        args=["Cl_stroma", "Cl_total"],
    )
    m.add_derived(
        "total_Cl_2",
        _squared,
        args=["Cl_total"],
    )
    m.add_derived(
        "total_K_2",
        _squared,
        args=["K_total"],
    )
    m.add_derived(
        "KEA3_reg",
        _reg_kea,
        args=["pH", "ATP", "KEA3_pH_reg", "KEA3_ATP_treshold"],
    )

    m.add_derived(
        "dG_K_ions",
        _dg_k,
        args=["K_lumen", "K_stroma", "delta_psi", "RT", "F"],
    )

    m.add_derived(
        "Cl_driving_force",
        _cl_driving_force,
        args=["delta_psi", "Cl_lumen", "Cl_stroma", "RT", "F"],
    )

    m.add_derived(
        "Keq_NDH1",
        _keq_ndh1,
        args=["pmf", "E^0_Fd", "F", "E^0_PQ", "pH", "dG_pH", "RT"],
    )

    m.add_reaction(
        "KEA3",
        _v_kea,
        args=["K_lumen", "protons_lumen", "K_stroma", "k_KEA", "protons", "KEA3_reg"],
        stoichiometry={
            "K_stroma": -1,
            "pH_lumen": 1 / bH,  # type: ignore
            "pH": -1 / buffering,
        },
    )

    m.add_reaction(
        "voltage_K_channel",
        _v_voltage_k_channel,
        args=[
            "delta_psi",
            "K_lumen",
            "K_stroma",
            "dG_K_ions",
            "perm_K",
            "K_delta_psi_treshold",
        ],
        stoichiometry={
            "K_stroma": 1,
            "delta_psi": Derived(
                fn=_neg_one_div,
                args=["volts_per_charge"],
            ),
        },
    )

    m.add_reaction(
        "VCCN1",
        _v_vccn1,
        args=[
            "Cl_stroma",
            "Cl_lumen",
            "Cl_driving_force",
            "delta_psi",
            "k_VCCN1",
            "VCCN_delta_psi_treshold",
        ],
        stoichiometry={
            "Cl_stroma": -1,
            "delta_psi": Derived(
                fn=_neg_one_div,
                args=["volts_per_charge"],
            ),
        },
    )

    m.add_reaction(
        "ClCe",
        _v_cl_ce,
        args=[
            "Cl_lumen",
            "Cl_stroma",
            "k_ClCe",
            n.pq_ox(),
            "Cl_driving_force",
            "ClCe_PQ",
        ],
        stoichiometry={
            "Cl_stroma": -1,
            "delta_psi": Derived(
                fn=_neg_one_div,
                args=["volts_per_charge"],
            ),
        },
    )

    m.add_reaction(
        "Cl_leak",
        _v_cl_leak,
        args=[
            "k_Cl_leak",
            "Cl_lumen",
            "Cl_stroma",
            n.pq_ox(),
            "Cl_leak_PQ",
            "total_Cl_2",
        ],
        stoichiometry={
            "Cl_stroma": 1,
            "delta_psi": Derived(
                fn=_one_div,
                args=["volts_per_charge"],
            ),
        },
    )

    m.add_reaction(
        "NDH1",
        _v_ndh1,
        args=["P700+FA-", n.fd_red(), n.pq_ox(), "pH_lumen", "k_NDH1"],
        stoichiometry={
            n.fd_ox(): 2,
            n.pq_ox(): -1,
            "pH_lumen": -4 / bH,  # type: ignore
            "pH": 4 / buffering,
            "delta_psi": Derived(
                fn=_four_div,
                args=["volts_per_charge"],
            ),
        },
    )

    if ClCe == "bi":
        m.remove_reaction("ClCe")
        m.add_parameter(
            "ClCe_ATP_threshold",
            0.2,
        )
        m.add_derived(
            "ClCe_activation",
            _cl_ce_activation,
            args=["ATP", "ClCe_ATP_threshold"],
        )

        m.add_reaction(
            "ClCe_bi",
            _cl_ce_bi,
            args=[
                "Cl_lumen",
                "Cl_stroma",
                "k_ClCe",
                "ClCe_activation",
            ],
            stoichiometry={"Cl_stroma": -1},
        )
        m.update_parameters({"k_ClCe": 0.5})

    if ClCe == "antiport":
        m.remove_reaction("ClCe")
        m.add_parameter(
            "ClCe_ATP_threshold",
            0.2,
        )
        m.add_derived(
            "ClCe_activation",
            _cl_ce_activation,
            args=["ATP", "ClCe_ATP_threshold"],
        )

        m.add_reaction(
            "ClCe_CH",
            _cl_ce_ch,
            args=[
                "Cl_lumen",
                "Cl_stroma",
                "k_ClCe",
                "ClCe_activation",
                "protons",
                "protons_lumen",
            ],
            stoichiometry={
                "Cl_stroma": 1,
                "pH_lumen": -1 / bH,  # type: ignore
                "pH": 1 / buffering,
            },
        )

    m.update_variables(
        {
            "pH": 7.5,
            "pH_lumen": 6.8,
        }
    )

    m.update_parameters(
        {
            "ClCe_PQ": 15.87880046767565,
            "Cl_leak_PQ": 14.92901445507139,
            "K_delta_psi_treshold": 0.081468073076241584,
            "VCCN_delta_psi_treshold": 0.080009009793326773,
            "k_Cl_leak": 25,
            "k_NDH1": 7.447430768265866,
            "perm_K": 1.6113135416150155,
            "k_VCCN1": 0.5,
            "k_ClCe": 0.5,
            "Cl_total": 50,
            "KEA3_pH_reg": 7.69,
            "KEA3_ATP_treshold": 0.26274793681796166,
            "k_KEA": 90,
            "pKreg": 7,
            "kActATPase": 0.001,
            "gamma0": 0.06260060801266355,
            "gamma1": 0.4053583123566203,
            "gamma2": 0.7040758738825374,
            "gamma3": 0.07834807781016208,
            "kH0": 5e8,
            "kH_Qslope": 3e10,
            "ksat_violaxanthin_deepoxidase": 6.1935954078503974,
            "kf_violaxanthin_deepoxidase": 0.0006091912188339879,
            "kf_zeaxanthin_epoxidase": 0.000106261953934132,
            "ksat_lhc_protonation": 6.2539066418842256,
            "kf_lhc_deprotonation": 0.015892570403695704,
            "kf_lhc_protonation": 0.15837051384170664,
            "pK0E": 5.960025833706074,
            "b": 1.8688304401249532,
            "K_total": 60,
            n.kh(n.violaxanthin_deepoxidase()): 4,
            n.kh(n.lhc_protonation()): 10,
            n.co2(): 0.013226,
            n.kf(n.ferredoxin_thioredoxin_reductase()): 0.8,
        }
    )

    return m


def get_ebeling_2026() -> Model:
    return _build_full_model(
        get_saadat2021(),
        ClCe="bi",
    )

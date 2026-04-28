from typing import Literal

import numpy as np
from mxlpy import Derived, Model
from mxlpy.surrogates import qss


def _ph(h: float) -> float:
    """Convert lumenal proton concentration to pH (conversion factor 2.5e-4)."""
    return -np.log10(h * 2.5e-4)


def _ph_inv(ph: float) -> float:
    """Convert pH to lumenal proton concentration (inverse of _ph)."""
    return 3.2e4 * 10**-ph


def _keq_qapq(
    f: float,
    e0_qa: float,
    e0_pq: float,
    p_h_st: float,
    r: float,
    t: float,
) -> float:
    """Equilibrium constant for QA -> PQ electron transfer, including stroma pH correction."""
    RT = r * t
    DG1 = -f * e0_qa
    DG2 = -2 * f * e0_pq + 2 * p_h_st * np.log(10) * RT
    DG0 = -2 * DG1 + DG2
    return np.exp(-DG0 / RT)


def _keq_cytb6f(
    p_h_lu: float,
    f: float,
    e0_pq: float,
    r: float,
    t: float,
    e0_pc: float,
    p_h_st: float,
) -> float:
    """Equilibriu constant of Cytochrome b6f."""
    RT = r * t
    DG1 = -2 * f * e0_pq + 2 * RT * np.log(10) * p_h_lu
    DG2 = -f * e0_pc
    DG3 = RT * np.log(10) * (p_h_st - p_h_lu)
    DG = -DG1 + 2 * DG2 + 2 * DG3
    return np.exp(-DG / RT)


def _keq_atp_synth(
    p_h_lu: float,
    dg_atp: float,
    p_h_st: float,
    r: float,
    t: float,
    pi: float,
) -> float:
    """Return equilibrium constant of ATP synthase.

    See Matuszynska et al 2016 or Ebenhöh et al. 2011,2014.
    """
    RT = r * t
    DG = dg_atp - np.log(10) * (p_h_st - p_h_lu) * (14 / 3) * RT
    return pi * np.exp(-DG / RT)


def _moiety_1s(
    x: float,
    x_total: float,
) -> float:
    """Calculate conservation relationship for one substrate.

    Used for creating derived variables that represent moiety conservation,
    such as calculating the free form of a species when you know the total.

    Parameters
    ----------
    x
        Concentration of one form of the species
    x_total
        Total concentration of all forms

    Returns
    -------
    Float
        Concentration of the other form (x_total - x)

    Examples
    --------
    >>> moiety_1s(0.3, 1.0)  # If total is 1.0 and one form is 0.3, other is 0.7
    0.7
    >>> # Example: If ATP + ADP = total_adenosine
    >>> moiety_1s(0.8, 1.5)  # [ADP] = total_adenosine - [ATP]
    0.7

    """
    return x_total - x


def _quencher(
    psb_s: float,
    zx: float,
    psb_sp: float,
    k_z_sat: float,
    gamma_0: float,
    gamma_1: float,
    gamma_2: float,
    gamma_3: float,
) -> float:
    """Quencher mechanism.

    accepts:
    psbS: fraction of non-protonated PsbS protein
    Vx: fraction of Violaxanthin
    """
    Zs = zx / (zx + k_z_sat)

    return (
        gamma_0 * (1 - Zs) * psb_s
        + gamma_1 * (1 - Zs) * psb_sp
        + gamma_2 * Zs * psb_sp
        + gamma_3 * Zs * psb_s
    )


def _fluorescence(
    q: float,
    b0: float,
    b2: float,
    k_h: float,
    k_f: float,
    k_p: float,
) -> float:
    """Fluorescence function."""
    return k_f / (k_h * q + k_f + k_p) * b0 + k_f / (k_h * q + k_f) * b2


def _psii(
    b1: float,
    k_p: float,
) -> float:
    """Reduction of PQ due to ps2."""
    return k_p * 0.5 * b1


def _two_divided_value(x: float) -> float:
    """Return 2/x; used for scaling 2-proton stoichiometry by buffering capacity."""
    return 2 / x


def _ptox(
    pqh_2: float,
    pfd: float,
    k_cytb6f: float,
    k_ptox: float,
    o2_ex: float,
    pq_tot: float,
    keq_cytb6f: float,
) -> float:
    """Oxidation of the PQ pool through cytochrome and PTOX."""
    kPFD = k_cytb6f * pfd
    k_ptox = k_ptox * o2_ex
    a1 = kPFD * keq_cytb6f / (keq_cytb6f + 1)
    a2 = kPFD / (keq_cytb6f + 1)
    return (a1 + k_ptox) * pqh_2 - a2 * (pq_tot - pqh_2)


def _four_divided_value(x: float) -> float:
    """Return 4/x; used for scaling 4-proton stoichiometry by buffering capacity."""
    return 4 / x


def _atp_synthase(
    atp_st: float,
    at_pase_ac: float,
    kf_at_psynth: float,
    keq_at_psynth: float,
    ap_tot: float,
) -> float:
    """Production of ATP by ATPsynthase."""
    return at_pase_ac * kf_at_psynth * (ap_tot - atp_st - atp_st / keq_at_psynth)


def _neg_fourteenthirds_divided_value(x: float) -> float:
    """Return -(14/3)/x; used for HPR proton stoichiometry of ATP synthase scaled by buffering capacity."""
    return -(14 / 3) / x


def _atp_synthase_activase(
    at_pase_ac: float,
    pfd: float,
    k_act_at_pase: float,
    k_deact_at_pase: float,
) -> float:
    """Activation of ATPsynthase by light."""
    switch = pfd > 0.0
    if switch:
        return k_act_at_pase * (1 - at_pase_ac)
    return -k_deact_at_pase * at_pase_ac


def _proton_leak(
    h_lu: float,
    k_leak: float,
    h_st: float,
) -> float:
    """Transmembrane proton leak."""
    return k_leak * (h_lu - h_st)


def _neg_divided_value(x: float) -> float:
    """Return -1/x; used for negative proton leak stoichiometry scaled by buffering capacity."""
    return -1 / x


def _atp_consumption(
    atp_st: float,
    k_at_pconsum: float,
) -> float:
    """ATP consuming reaction."""
    return k_at_pconsum * atp_st


def _xantophyll_cycle(
    vx: float,
    h_lu: float,
    nhx: float,
    k_p_h_sat_inv: float,
    k_dv: float,
    k_ez: float,
    x_tot: float,
) -> float:
    """Xanthophyll cycle."""
    a = h_lu**nhx / (h_lu**nhx + k_p_h_sat_inv**nhx)
    return k_dv * a * vx - k_ez * (x_tot - vx)


def _psbs_protonation(
    psb_s: float,
    h_lu: float,
    nhl: float,
    k_p_h_sat_lhc_inv: float,
    k_prot: float,
    k_deprot: float,
    psb_s_tot: float,
) -> float:
    """Protonation of PsbS protein."""
    a = h_lu**nhl / (h_lu**nhl + k_p_h_sat_lhc_inv**nhl)
    return k_prot * a * psb_s - k_deprot * (psb_s_tot - psb_s)


def _ps2states_2016a_matrix(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> tuple[float, float, float, float]:
    """PSII state populations (Matuszynska 2016 NPQ, matrix solver variant)."""
    b0 = pfd + k_pqh2 * pq_red / keq_qapq
    b1 = kh * quencher + kf
    b2 = kh * quencher + kf + kp

    M = np.array(
        [
            [-b0, b1, k_pqh2 * pq_ox, 0],  # B0
            [pfd, -b2, 0, 0],  # B1
            [0, 0, pfd, -b1],  # B3
            [1, 1, 1, 1],
        ]
    )

    A = np.array([0, 0, 0, psii_tot])
    b0, b1, b2, B3 = np.linalg.solve(M, A)
    return b0, b1, b2, B3


def _ps2states_2016a_analytical(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> tuple[float, float, float, float]:
    """PSII state populations (Matuszynska 2016 NPQ, analytical closed-form variant)."""
    x0 = kf**2
    x1 = kf * kp
    x2 = kh * quencher
    x3 = kp * x2
    x4 = 2 * x2
    x5 = kf * x4
    x6 = kh**2 * quencher**2
    x7 = keq_qapq * kp
    x8 = k_pqh2 * pq_ox
    x9 = keq_qapq * x8
    x10 = k_pqh2 * pq_red
    x11 = kf * x10
    x12 = kp * x10
    x13 = pfd * x9
    x14 = x10 * x2
    x15 = pfd * x7
    x16 = (
        keq_qapq * pfd * x1
        + x0 * x10
        + x1 * x10
        + x10 * x3
        + x10 * x6
        + x11 * x4
        + x15 * x2
    )
    x17 = psii_tot / (
        kf * x13
        + pfd**2 * x7
        + pfd * x11
        + pfd * x12
        + pfd * x14
        + x0 * x9
        + x1 * x9
        + x13 * x2
        + x16
        + x2 * x7 * x8
        + x5 * x9
        + x6 * x9
    )
    x18 = pfd * x17
    _B0 = x17 * x9 * (x0 + x1 + x3 + x5 + x6)
    _B1 = x18 * x9 * (kf + x2)
    _B2 = x16 * x17
    _B3 = x18 * (x11 + x12 + x14 + x15)
    return _B0, _B1, _B2, _B3


def b0(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> float:
    """Analytical PSII B0 state (dark-closed, Matuszynska 2016): fraction with PQ oxidised."""
    return (
        k_pqh2
        * keq_qapq
        * pq_ox
        * psii_tot
        * (
            kf**2
            + 2 * kf * kh * quencher
            + kf * kp
            + kh**2 * quencher**2
            + kh * kp * quencher
        )
        / (
            k_pqh2 * keq_qapq * kf**2 * pq_ox
            + 2 * k_pqh2 * keq_qapq * kf * kh * pq_ox * quencher
            + k_pqh2 * keq_qapq * kf * kp * pq_ox
            + k_pqh2 * keq_qapq * kf * pfd * pq_ox
            + k_pqh2 * keq_qapq * kh**2 * pq_ox * quencher**2
            + k_pqh2 * keq_qapq * kh * kp * pq_ox * quencher
            + k_pqh2 * keq_qapq * kh * pfd * pq_ox * quencher
            + k_pqh2 * kf**2 * pq_red
            + 2 * k_pqh2 * kf * kh * pq_red * quencher
            + k_pqh2 * kf * kp * pq_red
            + k_pqh2 * kf * pfd * pq_red
            + k_pqh2 * kh**2 * pq_red * quencher**2
            + k_pqh2 * kh * kp * pq_red * quencher
            + k_pqh2 * kh * pfd * pq_red * quencher
            + k_pqh2 * kp * pfd * pq_red
            + keq_qapq * kf * kp * pfd
            + keq_qapq * kh * kp * pfd * quencher
            + keq_qapq * kp * pfd**2
        )
    )


def b1(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> float:
    """Analytical PSII B1 state (light-open, Matuszynska 2016): fraction absorbing light with PQ oxidised."""
    return (
        k_pqh2
        * keq_qapq
        * pfd
        * pq_ox
        * psii_tot
        * (kf + kh * quencher)
        / (
            k_pqh2 * keq_qapq * kf**2 * pq_ox
            + 2 * k_pqh2 * keq_qapq * kf * kh * pq_ox * quencher
            + k_pqh2 * keq_qapq * kf * kp * pq_ox
            + k_pqh2 * keq_qapq * kf * pfd * pq_ox
            + k_pqh2 * keq_qapq * kh**2 * pq_ox * quencher**2
            + k_pqh2 * keq_qapq * kh * kp * pq_ox * quencher
            + k_pqh2 * keq_qapq * kh * pfd * pq_ox * quencher
            + k_pqh2 * kf**2 * pq_red
            + 2 * k_pqh2 * kf * kh * pq_red * quencher
            + k_pqh2 * kf * kp * pq_red
            + k_pqh2 * kf * pfd * pq_red
            + k_pqh2 * kh**2 * pq_red * quencher**2
            + k_pqh2 * kh * kp * pq_red * quencher
            + k_pqh2 * kh * pfd * pq_red * quencher
            + k_pqh2 * kp * pfd * pq_red
            + keq_qapq * kf * kp * pfd
            + keq_qapq * kh * kp * pfd * quencher
            + keq_qapq * kp * pfd**2
        )
    )


def b2(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> float:
    """Analytical PSII B2 state (dark-closed, Matuszynska 2016): fraction with PQ reduced."""
    return (
        psii_tot
        * (
            k_pqh2 * kf**2 * pq_red
            + 2 * k_pqh2 * kf * kh * pq_red * quencher
            + k_pqh2 * kf * kp * pq_red
            + k_pqh2 * kh**2 * pq_red * quencher**2
            + k_pqh2 * kh * kp * pq_red * quencher
            + keq_qapq * kf * kp * pfd
            + keq_qapq * kh * kp * pfd * quencher
        )
        / (
            k_pqh2 * keq_qapq * kf**2 * pq_ox
            + 2 * k_pqh2 * keq_qapq * kf * kh * pq_ox * quencher
            + k_pqh2 * keq_qapq * kf * kp * pq_ox
            + k_pqh2 * keq_qapq * kf * pfd * pq_ox
            + k_pqh2 * keq_qapq * kh**2 * pq_ox * quencher**2
            + k_pqh2 * keq_qapq * kh * kp * pq_ox * quencher
            + k_pqh2 * keq_qapq * kh * pfd * pq_ox * quencher
            + k_pqh2 * kf**2 * pq_red
            + 2 * k_pqh2 * kf * kh * pq_red * quencher
            + k_pqh2 * kf * kp * pq_red
            + k_pqh2 * kf * pfd * pq_red
            + k_pqh2 * kh**2 * pq_red * quencher**2
            + k_pqh2 * kh * kp * pq_red * quencher
            + k_pqh2 * kh * pfd * pq_red * quencher
            + k_pqh2 * kp * pfd * pq_red
            + keq_qapq * kf * kp * pfd
            + keq_qapq * kh * kp * pfd * quencher
            + keq_qapq * kp * pfd**2
        )
    )


def b3(
    pq_ox: float,
    pq_red: float,
    quencher: float,
    pfd: float,
    k_pqh2: float,
    keq_qapq: float,
    kh: float,
    kf: float,
    kp: float,
    psii_tot: float,
) -> float:
    """Analytical PSII B3 state (light-open, Matuszynska 2016): fraction absorbing light with PQ reduced."""
    return (
        pfd
        * psii_tot
        * (
            k_pqh2 * kf * pq_red
            + k_pqh2 * kh * pq_red * quencher
            + k_pqh2 * kp * pq_red
            + keq_qapq * kp * pfd
        )
        / (
            k_pqh2 * keq_qapq * kf**2 * pq_ox
            + 2 * k_pqh2 * keq_qapq * kf * kh * pq_ox * quencher
            + k_pqh2 * keq_qapq * kf * kp * pq_ox
            + k_pqh2 * keq_qapq * kf * pfd * pq_ox
            + k_pqh2 * keq_qapq * kh**2 * pq_ox * quencher**2
            + k_pqh2 * keq_qapq * kh * kp * pq_ox * quencher
            + k_pqh2 * keq_qapq * kh * pfd * pq_ox * quencher
            + k_pqh2 * kf**2 * pq_red
            + 2 * k_pqh2 * kf * kh * pq_red * quencher
            + k_pqh2 * kf * kp * pq_red
            + k_pqh2 * kf * pfd * pq_red
            + k_pqh2 * kh**2 * pq_red * quencher**2
            + k_pqh2 * kh * kp * pq_red * quencher
            + k_pqh2 * kh * pfd * pq_red * quencher
            + k_pqh2 * kp * pfd * pq_red
            + keq_qapq * kf * kp * pfd
            + keq_qapq * kh * kp * pfd * quencher
            + keq_qapq * kp * pfd**2
        )
    )


def get_matuszynska2016_npq(
    mode: Literal["matrix", "analytical", "analytical-split"] = "analytical",
) -> Model:
    """Build the Matuszynska 2016 NPQ model with PSII state computation in the chosen mode."""
    m = Model()

    m.add_variable("pq_red", initial_value=0)
    m.add_variable("protons", initial_value=6.32975752e-05)
    m.add_variable("vmax_atp_synthase", initial_value=0)
    m.add_variable("atp", initial_value=25.0)
    m.add_variable("psbs_de", initial_value=1)
    m.add_variable("vx", initial_value=1)
    m.add_parameter("PSII_tot", value=2.5)
    m.add_parameter("PQ_tot", value=20)
    m.add_parameter("AP_tot", value=50)
    m.add_parameter("PsbS_tot", value=1)
    m.add_parameter("X_tot", value=1)
    m.add_parameter("O2_ex", value=8)
    m.add_parameter("Pi", value=0.01)
    m.add_parameter("k_Cytb6f", value=0.104)
    m.add_parameter("k_ActATPase", value=0.01)
    m.add_parameter("k_DeactATPase", value=0.002)
    m.add_parameter("k_ATPsynth", value=20.0)
    m.add_parameter("k_ATPconsum", value=10.0)
    m.add_parameter("k_PQH2", value=250.0)
    m.add_parameter("k_H", value=5000000000.0)
    m.add_parameter("k_F", value=625000000.0)
    m.add_parameter("k_P", value=5000000000.0)
    m.add_parameter("k_PTOX", value=0.01)
    m.add_parameter("pH_st", value=7.8)
    m.add_parameter("k_leak", value=1000)
    m.add_parameter("b_H", value=100)
    m.add_parameter("hpr", value=4.666666666666667)
    m.add_parameter("k_DV", value=0.0024)
    m.add_parameter("k_EZ", value=0.00024)
    m.add_parameter("K_pHSat", value=5.8)
    m.add_parameter("nhx", value=5.0)
    m.add_parameter("K_ZSat", value=0.12)
    m.add_parameter("nhl", value=3)
    m.add_parameter("k_deprot", value=0.0096)
    m.add_parameter("k_prot", value=0.0096)
    m.add_parameter("K_pHSatLHC", value=5.8)
    m.add_parameter("gamma_0", value=0.1)
    m.add_parameter("gamma_1", value=0.25)
    m.add_parameter("gamma_2", value=0.6)
    m.add_parameter("gamma_3", value=0.15)
    m.add_parameter("F", value=96.485)
    m.add_parameter("R", value=0.0083)
    m.add_parameter("T", value=298)
    m.add_parameter("E0_QA", value=-0.14)
    m.add_parameter("E0_PQ", value=0.354)
    m.add_parameter("E0_PC", value=0.38)
    m.add_parameter("DG_ATP", value=30.6)
    m.add_parameter("PPFD", value=100)
    m.add_derived(
        "pH_lu",
        fn=_ph,
        args=["protons"],
    )
    m.add_derived(
        "H_st",
        fn=_ph_inv,
        args=["pH_st"],
    )
    m.add_derived(
        "K_pHSat_inv",
        fn=_ph_inv,
        args=["K_pHSat"],
    )
    m.add_derived(
        "K_pHSatLHC_inv",
        fn=_ph_inv,
        args=["K_pHSatLHC"],
    )
    m.add_derived(
        "K_QAPQ",
        fn=_keq_qapq,
        args=["F", "E0_QA", "E0_PQ", "pH_st", "R", "T"],
    )
    m.add_derived(
        "K_cytb6f",
        fn=_keq_cytb6f,
        args=["pH_lu", "F", "E0_PQ", "R", "T", "E0_PC", "pH_st"],
    )
    m.add_derived(
        "K_ATPsynth",
        fn=_keq_atp_synth,
        args=["pH_lu", "DG_ATP", "pH_st", "R", "T", "Pi"],
    )
    m.add_derived(
        "pq_ox",
        fn=_moiety_1s,
        args=["pq_red", "PQ_tot"],
    )
    m.add_derived(
        "adp",
        fn=_moiety_1s,
        args=["atp", "AP_tot"],
    )
    m.add_derived(
        "psbs_pr",
        fn=_moiety_1s,
        args=["psbs_de", "PsbS_tot"],
    )
    m.add_derived(
        "zx",
        fn=_moiety_1s,
        args=["vx", "X_tot"],
    )
    m.add_derived(
        "Q",
        fn=_quencher,
        args=[
            "psbs_de",
            "zx",
            "psbs_pr",
            "K_ZSat",
            "gamma_0",
            "gamma_1",
            "gamma_2",
            "gamma_3",
        ],
    )
    m.add_derived(
        "Fluo",
        fn=_fluorescence,
        args=["Q", "B0", "B2", "k_H", "k_F", "k_P"],
    )
    m.add_reaction(
        "v_PSII",
        fn=_psii,
        args=["B1", "k_P"],
        stoichiometry={
            "pq_red": 1,
            "protons": Derived(fn=_two_divided_value, args=["b_H"]),
        },
    )
    m.add_reaction(
        "v_PQ",
        fn=_ptox,
        args=[
            "pq_red",
            "PPFD",
            "k_Cytb6f",
            "k_PTOX",
            "O2_ex",
            "PQ_tot",
            "K_cytb6f",
        ],
        stoichiometry={
            "pq_red": -1,
            "protons": Derived(fn=_four_divided_value, args=["b_H"]),
        },
    )
    m.add_reaction(
        "atp_synthase",
        fn=_atp_synthase,
        args=["atp", "vmax_atp_synthase", "k_ATPsynth", "K_ATPsynth", "AP_tot"],
        stoichiometry={
            "atp": 1,
            "protons": Derived(fn=_neg_fourteenthirds_divided_value, args=["b_H"]),
        },
    )
    m.add_reaction(
        "atp_activase",
        fn=_atp_synthase_activase,
        args=["vmax_atp_synthase", "PPFD", "k_ActATPase", "k_DeactATPase"],
        stoichiometry={"vmax_atp_synthase": 1},
    )
    m.add_reaction(
        "proton_leak",
        fn=_proton_leak,
        args=["protons", "k_leak", "H_st"],
        stoichiometry={"protons": Derived(fn=_neg_divided_value, args=["b_H"])},
    )
    m.add_reaction(
        "v_ATPcons",
        fn=_atp_consumption,
        args=["atp", "k_ATPconsum"],
        stoichiometry={"atp": -1},
    )
    m.add_reaction(
        "v_Xcyc",
        fn=_xantophyll_cycle,
        args=["vx", "protons", "nhx", "K_pHSat_inv", "k_DV", "k_EZ", "X_tot"],
        stoichiometry={"vx": -1},
    )
    m.add_reaction(
        "v_PsbSP",
        fn=_psbs_protonation,
        args=[
            "psbs_de",
            "protons",
            "nhl",
            "K_pHSatLHC_inv",
            "k_prot",
            "k_deprot",
            "PsbS_tot",
        ],
        stoichiometry={"psbs_de": -1},
    )
    if mode == "matrix":
        m.add_surrogate(
            "ps2states",
            qss.Surrogate(
                model=_ps2states_2016a_matrix,
                args=[
                    "pq_ox",
                    "pq_red",
                    "Q",
                    "PPFD",
                    "k_PQH2",
                    "K_QAPQ",
                    "k_H",
                    "k_F",
                    "k_P",
                    "PSII_tot",
                ],
                outputs=["B0", "B1", "B2", "B3"],
            ),
        )
    elif mode == "analytical":
        m.add_surrogate(
            "ps2states",
            qss.Surrogate(
                model=_ps2states_2016a_analytical,
                args=[
                    "pq_ox",
                    "pq_red",
                    "Q",
                    "PPFD",
                    "k_PQH2",
                    "K_QAPQ",
                    "k_H",
                    "k_F",
                    "k_P",
                    "PSII_tot",
                ],
                outputs=["B0", "B1", "B2", "B3"],
            ),
        )
    else:
        m.add_derived(
            "B0",
            b0,
            args=[
                "pq_ox",
                "pq_red",
                "Q",
                "PPFD",
                "k_PQH2",
                "K_QAPQ",
                "k_H",
                "k_F",
                "k_P",
                "PSII_tot",
            ],
        )
        m.add_derived(
            "B1",
            b1,
            args=[
                "pq_ox",
                "pq_red",
                "Q",
                "PPFD",
                "k_PQH2",
                "K_QAPQ",
                "k_H",
                "k_F",
                "k_P",
                "PSII_tot",
            ],
        )
        m.add_derived(
            "B2",
            b2,
            args=[
                "pq_ox",
                "pq_red",
                "Q",
                "PPFD",
                "k_PQH2",
                "K_QAPQ",
                "k_H",
                "k_F",
                "k_P",
                "PSII_tot",
            ],
        )
        m.add_derived(
            "B3",
            b3,
            args=[
                "pq_ox",
                "pq_red",
                "Q",
                "PPFD",
                "k_PQH2",
                "K_QAPQ",
                "k_H",
                "k_F",
                "k_P",
                "PSII_tot",
            ],
        )
    return m

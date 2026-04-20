from typing import Literal

import numpy as np
from mxlpy import Derived, Model
from mxlpy.surrogates import qss


def _ph(h: float) -> float:
    return -np.log10(h * 2.5e-4)


def _ph_inv(ph: float) -> float:
    return 3.2e4 * 10**-ph


def _keq_qapq(
    F: float,
    E0_QA: float,
    E0_PQ: float,
    pH_st: float,
    R: float,
    T: float,
) -> float:
    RT = R * T
    DG1 = -F * E0_QA
    DG2 = -2 * F * E0_PQ + 2 * pH_st * np.log(10) * RT
    DG0 = -2 * DG1 + DG2
    return np.exp(-DG0 / RT)


def _keq_cytb6f(
    pH_lu: float,
    F: float,
    E0_PQ: float,
    R: float,
    T: float,
    E0_PC: float,
    pH_st: float,
) -> float:
    """Equilibriu constant of Cytochrome b6f."""
    RT = R * T
    DG1 = -2 * F * E0_PQ + 2 * RT * np.log(10) * pH_lu
    DG2 = -F * E0_PC
    DG3 = RT * np.log(10) * (pH_st - pH_lu)
    DG = -DG1 + 2 * DG2 + 2 * DG3
    return np.exp(-DG / RT)


def _keq_atp_synth(
    pH_lu: float,
    DG_ATP: float,
    pH_st: float,
    R: float,
    T: float,
    Pi: float,
) -> float:
    """Equilibrium constant of ATP synthase. For more
    information see Matuszynska et al 2016 or Ebenhöh et al. 2011,2014.
    """
    RT = R * T
    DG = DG_ATP - np.log(10) * (pH_st - pH_lu) * (14 / 3) * RT
    return Pi * np.exp(-DG / RT)


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
    psbS: float,
    Zx: float,
    PsbSP: float,
    K_ZSat: float,
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
    Zs = Zx / (Zx + K_ZSat)

    return (
        gamma_0 * (1 - Zs) * psbS
        + gamma_1 * (1 - Zs) * PsbSP
        + gamma_2 * Zs * PsbSP
        + gamma_3 * Zs * psbS
    )


def _fluorescence(
    Q: float,
    B0: float,
    B2: float,
    k_H: float,
    k_F: float,
    k_P: float,
) -> float:
    """Fluorescence function."""
    return k_F / (k_H * Q + k_F + k_P) * B0 + k_F / (k_H * Q + k_F) * B2


def _psii(
    B1: float,
    k_P: float,
) -> float:
    """Reduction of PQ due to ps2."""
    return k_P * 0.5 * B1


def _two_divided_value(x: float) -> float:
    return 2 / x


def _ptox(
    PQH_2: float,
    pfd: float,
    k_Cytb6f: float,
    k_PTOX: float,
    O2_ex: float,
    PQ_tot: float,
    K_cytb6f: float,
) -> float:
    """Oxidation of the PQ pool through cytochrome and PTOX."""
    kPFD = k_Cytb6f * pfd
    k_PTOX = k_PTOX * O2_ex
    a1 = kPFD * K_cytb6f / (K_cytb6f + 1)
    a2 = kPFD / (K_cytb6f + 1)
    return (a1 + k_PTOX) * PQH_2 - a2 * (PQ_tot - PQH_2)


def _four_divided_value(x: float) -> float:
    return 4 / x


def _atp_synthase(
    ATP_st: float,
    ATPase_ac: float,
    k_ATPsynth: float,
    K_ATPsynth: float,
    AP_tot: float,
) -> float:
    """Production of ATP by ATPsynthase."""
    return ATPase_ac * k_ATPsynth * (AP_tot - ATP_st - ATP_st / K_ATPsynth)


def _neg_fourteenthirds_divided_value(x: float) -> float:
    return -(14 / 3) / x


def _atp_synthase_activase(
    ATPase_ac: float,
    pfd: float,
    k_ActATPase: float,
    k_DeactATPase: float,
) -> float:
    """Activation of ATPsynthase by light."""
    switch = pfd > 0.0
    if switch:
        return k_ActATPase * (1 - ATPase_ac)
    return -k_DeactATPase * ATPase_ac


def _proton_leak(
    H_lu: float,
    k_leak: float,
    H_st: float,
) -> float:
    """Transmembrane proton leak."""
    return k_leak * (H_lu - H_st)


def _neg_divided_value(x: float) -> float:
    return -1 / x


def _atp_consumption(
    ATP_st: float,
    k_ATPconsum: float,
) -> float:
    """ATP consuming reaction."""
    return k_ATPconsum * ATP_st


def _xantophyll_cycle(
    Vx: float,
    H_lu: float,
    nhx: float,
    K_pHSat_inv: float,
    k_DV: float,
    k_EZ: float,
    X_tot: float,
) -> float:
    """Xanthophyll cycle."""
    a = H_lu**nhx / (H_lu**nhx + K_pHSat_inv**nhx)
    return k_DV * a * Vx - k_EZ * (X_tot - Vx)


def _psbs_protonation(
    psbS: float,
    H_lu: float,
    nhl: float,
    K_pHSatLHC_inv: float,
    k_prot: float,
    k_deprot: float,
    PsbS_tot: float,
) -> float:
    """Protonation of PsbS protein."""
    a = H_lu**nhl / (H_lu**nhl + K_pHSatLHC_inv**nhl)
    return k_prot * a * psbS - k_deprot * (PsbS_tot - psbS)


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
    B0, B1, B2, B3 = np.linalg.solve(M, A)
    return B0, B1, B2, B3


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

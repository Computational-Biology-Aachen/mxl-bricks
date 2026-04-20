"""name

EC FIXME

Equilibrator
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np
from mxlpy import Derived, Model
from mxlpy.surrogates import qss

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, mass_action_2s, value
from mxlbricks.utils import (
    default_name,
    default_par,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from collections.abc import Iterable


def _two_div_by(x: float) -> float:
    return 2.0 / x


def _keq_pcp700(
    e0_pc: float,
    f: float,
    eo_p700: float,
    rt: float,
) -> float:
    dg1 = -e0_pc * f
    dg2 = -eo_p700 * f
    dg = -dg1 + dg2
    return math.exp(-dg / rt)


def _keq_faf_d(
    e0_fa: float,
    f: float,
    e0_fd: float,
    rt: float,
) -> float:
    dg1 = -e0_fa * f
    dg2 = -e0_fd * f
    dg = -dg1 + dg2
    return math.exp(-dg / rt)


def _rate_ps1(
    a: float,
    ps2cs: float,
    pfd: float,
) -> float:
    return (1 - ps2cs) * pfd * a


def _rate_ps2(
    b1: float,
    k2: float,
) -> float:
    return 0.5 * k2 * b1


def _ps1states_2019(
    pc_px: float,
    pc_red: float,
    fd_ox: float,
    fd_red: float,
    ps2cs: float,
    psi_tot: float,
    k_fd_red: float,
    keq_fafd: float,
    keq_pcp700: float,
    k_pc_ox: float,
    pfd: float,
) -> float:
    """QSSA calculates open state of PSI
    depends on reduction states of plastocyanin and ferredoxin
    C = [PC], F = [Fd] (ox. forms)
    """
    L = (1 - ps2cs) * pfd
    return psi_tot / (
        1
        + L / (k_fd_red * fd_ox)
        + (1 + fd_red / (keq_fafd * fd_ox))
        * (pc_px / (keq_pcp700 * pc_red) + L / (k_pc_ox * pc_red))
    )


def _ps1states_2021_surrogate(
    pc_ox: float,
    pc_red: float,
    fd_ox: float,
    fd_red: float,
    ps2cs: float,
    ps1_tot: float,
    k_fd_red: float,
    keq_f: float,
    keq_c: float,
    k_pc_ox: float,
    pfd: float,
    k0: float,
    o2: float,
) -> tuple[float, float, float]:
    """QSSA calculates open state of PSI
    depends on reduction states of plastocyanin and ferredoxin
    C = [PC], F = [Fd] (ox. forms)
    """
    kLI = (1 - ps2cs) * pfd

    y0 = (
        keq_c
        * keq_f
        * pc_red
        * ps1_tot
        * k_pc_ox
        * (fd_ox * k_fd_red + o2 * k0)
        / (
            fd_ox * keq_c * keq_f * pc_red * k_fd_red * k_pc_ox
            + fd_ox * keq_f * k_fd_red * (keq_c * kLI + pc_ox * k_pc_ox)
            + fd_red * k_fd_red * (keq_c * kLI + pc_ox * k_pc_ox)
            + keq_c * keq_f * o2 * pc_red * k0 * k_pc_ox
            + keq_c * keq_f * pc_red * kLI * k_pc_ox
            + keq_f * o2 * k0 * (keq_c * kLI + pc_ox * k_pc_ox)
        )
    )

    y1 = (
        ps1_tot
        * (
            fd_red * k_fd_red * (keq_c * kLI + pc_ox * k_pc_ox)
            + keq_c * keq_f * pc_red * kLI * k_pc_ox
        )
        / (
            fd_ox * keq_c * keq_f * pc_red * k_fd_red * k_pc_ox
            + fd_ox * keq_f * k_fd_red * (keq_c * kLI + pc_ox * k_pc_ox)
            + fd_red * k_fd_red * (keq_c * kLI + pc_ox * k_pc_ox)
            + keq_c * keq_f * o2 * pc_red * k0 * k_pc_ox
            + keq_c * keq_f * pc_red * kLI * k_pc_ox
            + keq_f * o2 * k0 * (keq_c * kLI + pc_ox * k_pc_ox)
        )
    )
    y2 = ps1_tot - y0 - y1

    return y0, y1, y2


def _ps2_crosssection(
    lhc: float,
    static_ant_ii: float,
    static_ant_i: float,
) -> float:
    return static_ant_ii + (1 - static_ant_ii - static_ant_i) * lhc


def _ps2states_2016_npq_matrix(
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
) -> Iterable[float]:
    b0 = pfd + k_pqh2 * pq_red / keq_qapq
    b1 = kh * quencher + kf
    b2 = kh * quencher + kf + kp

    state_matrix = np.array(
        [
            [-b0, b1, k_pqh2 * pq_ox, 0],  # B0
            [pfd, -b2, 0, 0],  # B1
            [0, 0, pfd, -b1],  # B3
            [1, 1, 1, 1],
        ]
    )

    a = np.array([0, 0, 0, psii_tot])
    return np.linalg.solve(state_matrix, a)


def _ps2states_2016_phd_matrix(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    q: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> Iterable[float]:
    absorbed = ps2cs * pfd
    kH = k_h0 + _kh * q
    k3p = k_pq_red * pq_ox
    k3m = k_pq_red * pq_red / keq_pq_red

    state_matrix = np.array(
        [
            [-absorbed - k3m, kH + k_f, k3p, 0],
            [absorbed, -(kH + k_f + k2), 0, 0],
            [0, 0, absorbed, -(kH + k_f)],
            [1, 1, 1, 1],
        ],
        dtype=float,
    )
    a = np.array([0, 0, 0, psii_tot])

    return np.linalg.solve(state_matrix, a)


def _ps2states_2016_npq_surrogate(
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


def _ps2states_2016_phd_surrogate(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    quencher: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> tuple[float, float, float, float]:
    x0 = k_f**2
    x1 = k_h0**2
    x2 = k2 * k_f
    x3 = k2 * k_h0
    x4 = 2 * k_f
    x5 = k_h0 * x4
    x6 = _kh * quencher
    x7 = k2 * x6
    x8 = x4 * x6
    x9 = 2 * x6
    x10 = k_h0 * x9
    x11 = _kh**2 * quencher**2
    x12 = k2 * keq_pq_red
    x13 = k_pq_red * keq_pq_red * pq_ox
    x14 = k_pq_red * pq_red
    x15 = k2 * x14
    x16 = pfd * ps2cs
    x17 = k_f * x14
    x18 = k_h0 * x14
    x19 = x14 * x6
    x20 = x13 * x16
    x21 = keq_pq_red * x16
    x22 = (
        x0 * x14
        + x1 * x14
        + x11 * x14
        + x14 * x2
        + x14 * x3
        + x14 * x5
        + x14 * x7
        + x14 * x8
        + x18 * x9
        + x2 * x21
        + x21 * x3
        + x21 * x7
    )
    x23 = psii_tot / (
        k_f * x20
        + k_h0 * x20
        + pfd**2 * ps2cs**2 * x12
        + x0 * x13
        + x1 * x13
        + x10 * x13
        + x11 * x13
        + x13 * x2
        + x13 * x3
        + x13 * x5
        + x13 * x7
        + x13 * x8
        + x15 * x16
        + x16 * x17
        + x16 * x18
        + x16 * x19
        + x20 * x6
        + x22
    )
    x24 = x16 * x23
    _B0 = x13 * x23 * (x0 + x1 + x10 + x11 + x2 + x3 + x5 + x7 + x8)
    _B1 = x13 * x24 * (k_f + k_h0 + x6)
    _B2 = x22 * x23
    _B3 = x24 * (x12 * x16 + x15 + x17 + x18 + x19)
    return _B0, _B1, _B2, _B3


def _b0_npq(
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


def _b1_npq(
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


def _b2_npq(
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


def _b3_npq(
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


def _b0_phd(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    quencher: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> float:
    return (
        k_pq_red
        * keq_pq_red
        * pq_ox
        * psii_tot
        * (
            _kh**2 * quencher**2
            + _kh * k2 * quencher
            + 2 * _kh * k_f * quencher
            + 2 * _kh * k_h0 * quencher
            + k2 * k_f
            + k2 * k_h0
            + k_f**2
            + 2 * k_f * k_h0
            + k_h0**2
        )
        / (
            _kh**2 * k_pq_red * keq_pq_red * pq_ox * quencher**2
            + _kh**2 * k_pq_red * pq_red * quencher**2
            + _kh * k2 * k_pq_red * keq_pq_red * pq_ox * quencher
            + _kh * k2 * k_pq_red * pq_red * quencher
            + _kh * k2 * keq_pq_red * pfd * ps2cs * quencher
            + 2 * _kh * k_f * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_f * k_pq_red * pq_red * quencher
            + 2 * _kh * k_h0 * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_h0 * k_pq_red * pq_red * quencher
            + _kh * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs * quencher
            + _kh * k_pq_red * pfd * pq_red * ps2cs * quencher
            + k2 * k_f * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_f * k_pq_red * pq_red
            + k2 * k_f * keq_pq_red * pfd * ps2cs
            + k2 * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_h0 * k_pq_red * pq_red
            + k2 * k_h0 * keq_pq_red * pfd * ps2cs
            + k2 * k_pq_red * pfd * pq_red * ps2cs
            + k2 * keq_pq_red * pfd**2 * ps2cs**2
            + k_f**2 * k_pq_red * keq_pq_red * pq_ox
            + k_f**2 * k_pq_red * pq_red
            + 2 * k_f * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + 2 * k_f * k_h0 * k_pq_red * pq_red
            + k_f * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_f * k_pq_red * pfd * pq_red * ps2cs
            + k_h0**2 * k_pq_red * keq_pq_red * pq_ox
            + k_h0**2 * k_pq_red * pq_red
            + k_h0 * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_h0 * k_pq_red * pfd * pq_red * ps2cs
        )
    )


def _b1_phd(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    quencher: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> float:
    return (
        k_pq_red
        * keq_pq_red
        * pfd
        * pq_ox
        * ps2cs
        * psii_tot
        * (_kh * quencher + k_f + k_h0)
        / (
            _kh**2 * k_pq_red * keq_pq_red * pq_ox * quencher**2
            + _kh**2 * k_pq_red * pq_red * quencher**2
            + _kh * k2 * k_pq_red * keq_pq_red * pq_ox * quencher
            + _kh * k2 * k_pq_red * pq_red * quencher
            + _kh * k2 * keq_pq_red * pfd * ps2cs * quencher
            + 2 * _kh * k_f * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_f * k_pq_red * pq_red * quencher
            + 2 * _kh * k_h0 * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_h0 * k_pq_red * pq_red * quencher
            + _kh * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs * quencher
            + _kh * k_pq_red * pfd * pq_red * ps2cs * quencher
            + k2 * k_f * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_f * k_pq_red * pq_red
            + k2 * k_f * keq_pq_red * pfd * ps2cs
            + k2 * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_h0 * k_pq_red * pq_red
            + k2 * k_h0 * keq_pq_red * pfd * ps2cs
            + k2 * k_pq_red * pfd * pq_red * ps2cs
            + k2 * keq_pq_red * pfd**2 * ps2cs**2
            + k_f**2 * k_pq_red * keq_pq_red * pq_ox
            + k_f**2 * k_pq_red * pq_red
            + 2 * k_f * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + 2 * k_f * k_h0 * k_pq_red * pq_red
            + k_f * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_f * k_pq_red * pfd * pq_red * ps2cs
            + k_h0**2 * k_pq_red * keq_pq_red * pq_ox
            + k_h0**2 * k_pq_red * pq_red
            + k_h0 * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_h0 * k_pq_red * pfd * pq_red * ps2cs
        )
    )


def _b2_phd(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    quencher: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> float:
    return (
        psii_tot
        * (
            _kh**2 * k_pq_red * pq_red * quencher**2
            + _kh * k2 * k_pq_red * pq_red * quencher
            + _kh * k2 * keq_pq_red * pfd * ps2cs * quencher
            + 2 * _kh * k_f * k_pq_red * pq_red * quencher
            + 2 * _kh * k_h0 * k_pq_red * pq_red * quencher
            + k2 * k_f * k_pq_red * pq_red
            + k2 * k_f * keq_pq_red * pfd * ps2cs
            + k2 * k_h0 * k_pq_red * pq_red
            + k2 * k_h0 * keq_pq_red * pfd * ps2cs
            + k_f**2 * k_pq_red * pq_red
            + 2 * k_f * k_h0 * k_pq_red * pq_red
            + k_h0**2 * k_pq_red * pq_red
        )
        / (
            _kh**2 * k_pq_red * keq_pq_red * pq_ox * quencher**2
            + _kh**2 * k_pq_red * pq_red * quencher**2
            + _kh * k2 * k_pq_red * keq_pq_red * pq_ox * quencher
            + _kh * k2 * k_pq_red * pq_red * quencher
            + _kh * k2 * keq_pq_red * pfd * ps2cs * quencher
            + 2 * _kh * k_f * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_f * k_pq_red * pq_red * quencher
            + 2 * _kh * k_h0 * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_h0 * k_pq_red * pq_red * quencher
            + _kh * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs * quencher
            + _kh * k_pq_red * pfd * pq_red * ps2cs * quencher
            + k2 * k_f * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_f * k_pq_red * pq_red
            + k2 * k_f * keq_pq_red * pfd * ps2cs
            + k2 * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_h0 * k_pq_red * pq_red
            + k2 * k_h0 * keq_pq_red * pfd * ps2cs
            + k2 * k_pq_red * pfd * pq_red * ps2cs
            + k2 * keq_pq_red * pfd**2 * ps2cs**2
            + k_f**2 * k_pq_red * keq_pq_red * pq_ox
            + k_f**2 * k_pq_red * pq_red
            + 2 * k_f * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + 2 * k_f * k_h0 * k_pq_red * pq_red
            + k_f * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_f * k_pq_red * pfd * pq_red * ps2cs
            + k_h0**2 * k_pq_red * keq_pq_red * pq_ox
            + k_h0**2 * k_pq_red * pq_red
            + k_h0 * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_h0 * k_pq_red * pfd * pq_red * ps2cs
        )
    )


def _b3_phd(
    pq_ox: float,
    pq_red: float,
    ps2cs: float,
    quencher: float,
    psii_tot: float,
    k2: float,
    k_f: float,
    _kh: float,
    keq_pq_red: float,
    k_pq_red: float,
    pfd: float,
    k_h0: float,
) -> float:
    return (
        pfd
        * ps2cs
        * psii_tot
        * (
            _kh * k_pq_red * pq_red * quencher
            + k2 * k_pq_red * pq_red
            + k2 * keq_pq_red * pfd * ps2cs
            + k_f * k_pq_red * pq_red
            + k_h0 * k_pq_red * pq_red
        )
        / (
            _kh**2 * k_pq_red * keq_pq_red * pq_ox * quencher**2
            + _kh**2 * k_pq_red * pq_red * quencher**2
            + _kh * k2 * k_pq_red * keq_pq_red * pq_ox * quencher
            + _kh * k2 * k_pq_red * pq_red * quencher
            + _kh * k2 * keq_pq_red * pfd * ps2cs * quencher
            + 2 * _kh * k_f * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_f * k_pq_red * pq_red * quencher
            + 2 * _kh * k_h0 * k_pq_red * keq_pq_red * pq_ox * quencher
            + 2 * _kh * k_h0 * k_pq_red * pq_red * quencher
            + _kh * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs * quencher
            + _kh * k_pq_red * pfd * pq_red * ps2cs * quencher
            + k2 * k_f * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_f * k_pq_red * pq_red
            + k2 * k_f * keq_pq_red * pfd * ps2cs
            + k2 * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + k2 * k_h0 * k_pq_red * pq_red
            + k2 * k_h0 * keq_pq_red * pfd * ps2cs
            + k2 * k_pq_red * pfd * pq_red * ps2cs
            + k2 * keq_pq_red * pfd**2 * ps2cs**2
            + k_f**2 * k_pq_red * keq_pq_red * pq_ox
            + k_f**2 * k_pq_red * pq_red
            + 2 * k_f * k_h0 * k_pq_red * keq_pq_red * pq_ox
            + 2 * k_f * k_h0 * k_pq_red * pq_red
            + k_f * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_f * k_pq_red * pfd * pq_red * ps2cs
            + k_h0**2 * k_pq_red * keq_pq_red * pq_ox
            + k_h0**2 * k_pq_red * pq_red
            + k_h0 * k_pq_red * keq_pq_red * pfd * pq_ox * ps2cs
            + k_h0 * k_pq_red * pfd * pq_red * ps2cs
        )
    )


def add_ps2_cross_section(
    model: Model,
    lhc: str | None = None,
    static_ant_i: str | None = None,
    static_ant_ii: str | None = None,
) -> Model:
    model.add_derived(
        name=n.ps2cs(),
        fn=_ps2_crosssection,
        args=[
            default_name(lhc, n.lhc),
            default_par(model, par=static_ant_ii, name="staticAntII", value=0.1),
            default_par(model, par=static_ant_i, name="staticAntI", value=0.37),
        ],
    )
    return model


def add_psii(
    model: Model,
    mode: Literal["matrix", "analytical", "analytical-split"] = "analytical",
    *,
    rxn: str,
    pq_ox: str,
    pq_red: str,
    ps2cs: str,
    quencher: str,
    pfd: str,
    h_lumen: str,
    b0: str,
    b1: str,
    b2: str,
    b3: str,
) -> Model:
    args = [
        pq_ox,
        pq_red,
        ps2cs,
        quencher,
        "PSII_total",
        "k2",
        "kF",
        "kH",
        n.keq(pq_red),
        "kPQred",
        pfd,
        "kH0",
    ]

    if mode == "matrix":
        model.add_surrogate(
            "ps2states",
            surrogate=qss.Surrogate(
                model=_ps2states_2016_phd_matrix,
                args=args,
                outputs=[
                    b0,
                    b1,
                    b2,
                    b3,
                ],
            ),
        )
    elif mode == "analytical":
        model.add_surrogate(
            "ps2states",
            surrogate=qss.Surrogate(
                model=_ps2states_2016_phd_surrogate,
                args=args,
                outputs=[
                    b0,
                    b1,
                    b2,
                    b3,
                ],
            ),
        )
    elif mode == "analytical-split":
        model.add_derived("B0", _b0_phd, args=args)
        model.add_derived("B1", _b1_phd, args=args)
        model.add_derived("B2", _b2_phd, args=args)
        model.add_derived("B3", _b3_phd, args=args)
    else:
        msg = (
            f"Unknown mode {mode}. Allowed modes are\n"
            "- matrix\n"
            "- analytical\n"
            "- analytical-split\n"
        )
        raise KeyError(msg)

    model.add_reaction(
        name=rxn,
        fn=_rate_ps2,
        stoichiometry={
            pq_ox: -1,
            h_lumen: Derived(fn=_two_div_by, args=["bH"]),
        },
        args=[
            b1,
            "k2",
        ],
    )

    return model


def add_psi_2019(
    model: Model,
    *,
    rxn: str,
    ps2cs: str,
    pfd: str,
    pc_ox: str,
    pc_red: str,
    fd_ox: str,
    fd_red: str,
    a1: str,
    keq_pcp700: str,
    keq_fd_red: str,
) -> Model:
    model.add_derived(
        name=n.a1(),
        fn=_ps1states_2019,
        args=[
            pc_ox,
            pc_red,
            fd_ox,
            fd_red,
            ps2cs,
            "PSI_total",
            "kFdred",
            keq_fd_red,
            keq_pcp700,
            "kPCox",
            pfd,
        ],
    )
    model.add_reaction(
        name=rxn,
        fn=_rate_ps1,
        stoichiometry={
            fd_ox: -1,
            pc_ox: 1,
        },
        args=[
            a1,
            ps2cs,
            pfd,
        ],
    )
    return model


def add_psi_2021(
    model: Model,
    *,
    rxn: str,
    ps2cs: str,
    pfd: str,
    o2_lumen: str,
    pc_ox: str,
    pc_red: str,
    fd_ox: str,
    fd_red: str,
    a0: str,
    a1: str,
    a2: str,
    keq_pcp700: str,
    keq_fd_red: str,
    k_mehler: str,
) -> Model:
    model.add_surrogate(
        "ps1states",
        surrogate=qss.Surrogate(
            model=_ps1states_2021_surrogate,
            args=[
                pc_ox,
                pc_red,
                fd_ox,
                fd_red,
                ps2cs,
                "PSI_total",
                "kFdred",
                keq_fd_red,
                keq_pcp700,
                "kPCox",
                pfd,
                k_mehler,
                o2_lumen,
            ],
            outputs=[
                a0,
                a1,
                a2,
            ],
        ),
    )
    model.add_reaction(
        name=rxn,
        fn=_rate_ps1,
        stoichiometry={
            pc_ox: 1,
        },
        args=[
            a0,
            ps2cs,
            pfd,
        ],
    )
    return model


def add_mehler(
    model: Model,
    *,
    rxn: str,
    o2_lumen: str,
    h2o2: str,
    a1: str,
    convf: str,
    k_mehler: str,
) -> Model:
    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry={
            h2o2: Derived(fn=value, args=[convf]),
        },
        args=[
            a1,
            o2_lumen,
            k_mehler,
        ],
    )
    return model


def add_photosystems(
    model: Model,
    mode: Literal["matrix", "analytical", "analytical-split"] = "analytical",
    *,
    rxn_psii: str | None = None,
    rxn_psi: str | None = None,
    rxn_mehler: str | None = None,
    pq_ox: str | None = None,
    pq_red: str | None = None,
    ps2cs: str | None = None,
    quencher: str | None = None,
    pfd: str | None = None,
    o2_lumen: str | None = None,
    h_lumen: str | None = None,
    h2o2: str | None = None,
    pc_ox: str | None = None,
    pc_red: str | None = None,
    fd_ox: str | None = None,
    fd_red: str | None = None,
    a0: str | None = None,
    a1: str | None = None,
    a2: str | None = None,
    b0: str | None = None,
    b1: str | None = None,
    b2: str | None = None,
    b3: str | None = None,
    mehler: bool,
    convf: str | None = None,
) -> Model:
    """PSII: 2 H2O + 2 PQ + 4 H_stroma -> O2 + 2 PQH2 + 4 H_lumen
    PSI: Fd_ox + PC_red -> Fd_red + PC_ox
    """
    pq_ox = default_name(pq_ox, n.pq_ox)
    pq_red = default_name(pq_red, n.pq_red)
    ps2cs = default_name(ps2cs, n.ps2cs)
    quencher = default_name(quencher, n.quencher)
    pfd = default_name(pfd, n.pfd)
    o2_lumen = default_name(h_lumen, lambda: n.o2("_lumen"))
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    h2o2 = default_name(h2o2, n.h2o2)
    pc_ox = default_name(pc_ox, n.pc_ox)
    pc_red = default_name(pc_red, n.pc_red)
    fd_ox = default_name(fd_ox, n.fd_ox)
    fd_red = default_name(fd_red, n.fd_red)
    a0 = default_name(a0, n.a0)
    a1 = default_name(a1, n.a1)
    a2 = default_name(a2, n.a2)
    b0 = default_name(b0, n.b0)
    b1 = default_name(b1, n.b1)
    b2 = default_name(b2, n.b2)
    b3 = default_name(b3, n.b3)

    model.add_parameter("PSII_total", 2.5)
    model.add_parameter("PSI_total", 2.5)
    model.add_parameter("kH0", 500000000.0)
    model.add_parameter("kPQred", 250.0)
    model.add_parameter("kPCox", 2500.0)
    model.add_parameter("kFdred", 250000.0)
    model.add_parameter("k2", 5000000000.0)
    model.add_parameter("kH", 5000000000.0)
    model.add_parameter("kF", 625000000.0)
    convf = static(model, n.convf(), 3.2e-2) if convf is None else convf

    model.add_derived(
        keq_pcp700 := n.keq("PCP700"),
        _keq_pcp700,
        args=["E^0_PC", "F", "E^0_P700", "RT"],
    )
    model.add_derived(
        keq_fd_red := n.keq(n.ferredoxin_reductase()),
        _keq_faf_d,
        args=["E^0_FA", "F", "E^0_Fd", "RT"],
    )

    add_psii(
        model,
        mode=mode,
        rxn=default_name(rxn_psii, n.ps2),
        pq_ox=pq_ox,
        pq_red=pq_red,
        ps2cs=ps2cs,
        quencher=quencher,
        pfd=pfd,
        h_lumen=h_lumen,
        b0=b0,
        b1=b1,
        b2=b2,
        b3=b3,
    )

    if not mehler:
        add_psi_2019(
            model,
            rxn=default_name(rxn_psi, n.ps1),
            ps2cs=ps2cs,
            pfd=pfd,
            pc_ox=pc_ox,
            pc_red=pc_red,
            fd_ox=fd_ox,
            fd_red=fd_red,
            a1=a1,
            keq_pcp700=keq_pcp700,
            keq_fd_red=keq_fd_red,
        )
    else:
        model.add_parameter(k_mehler := "kMehler", 1.0)
        add_psi_2021(
            model,
            rxn=default_name(rxn_psi, n.ps1),
            ps2cs=ps2cs,
            pfd=pfd,
            o2_lumen=o2_lumen,
            pc_ox=pc_ox,
            pc_red=pc_red,
            fd_ox=fd_ox,
            fd_red=fd_red,
            a0=a0,
            a1=a1,
            a2=a2,
            keq_pcp700=keq_pcp700,
            keq_fd_red=keq_fd_red,
            k_mehler=k_mehler,
        )
        add_mehler(
            model,
            rxn=default_name(rxn_mehler, n.mehler),
            o2_lumen=o2_lumen,
            h2o2=h2o2,
            a1=a1,
            convf=convf,
            k_mehler=k_mehler,
        )
    return model


def add_energy_production(
    model: Model,
    *,
    rxn: str | None = None,
    energy: str | None = None,
    pfd: str | None = None,
) -> Model:
    rxn = default_name(rxn, n.ps2)
    energy = default_name(energy, n.energy)
    pfd = default_name(pfd, n.pfd)

    model.add_parameter(k := n.kcat(pfd), 1 / 145)  # Fitted
    model.add_parameter(pfd, 700)

    model.add_reaction(
        n.petc(),
        mass_action_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                # Substrates
                # Products
                energy: 1,
            },
        ),
        args=[
            pfd,
            k,
        ],
    )
    return model

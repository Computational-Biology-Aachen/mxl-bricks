"""Reusable rate-law and algebraic functions for mxlpy models."""

from __future__ import annotations

###############################################################################
# General-purpose functions
###############################################################################


def value(x: float) -> float:
    """Return x unchanged."""
    return x


def neg(x: float) -> float:
    """Return -x."""
    return -x


def minus(x: float, y: float) -> float:
    """Return x - y."""
    return x - y


def mul(x: float, y: float) -> float:
    """Return x * y."""
    return x * y


def div(x: float, y: float) -> float:
    """Return x / y."""
    return x / y


def one_div(x: float) -> float:
    """Return 1 / x."""
    return 1.0 / x


def neg_div(x: float, y: float) -> float:
    """Return -x / y."""
    return -x / y


def twice(x: float) -> float:
    """Return 2 * x."""
    return x * 2


def proportional(x: float, y: float) -> float:
    """Return x * y."""
    return x * y


def diffusion(inside: float, outside: float, k: float) -> float:
    """Diffusion rate proportional to concentration gradient."""
    return k * (outside - inside)


###############################################################################
# Common algebraic modules
###############################################################################


def moiety_1(concentration: float, total: float) -> float:
    """Conservation moiety: total - concentration."""
    return total - concentration


def moiety_2(x1: float, x2: float, total: float) -> float:
    """Conservation moiety: total - x1 - x2."""
    return total - x1 - x2


###############################################################################
# Common rate laws
###############################################################################


def mass_action_1s(s1: float, k_fwd: float) -> float:
    """Mass-action rate for one substrate."""
    return k_fwd * s1


def mass_action_2s(s1: float, s2: float, k_fwd: float) -> float:
    """Mass-action rate for two substrates."""
    return k_fwd * s1 * s2


def mass_action_3s(s1: float, s2: float, s3: float, k_fwd: float) -> float:
    """Mass-action rate for three substrates."""
    return k_fwd * s1 * s2 * s3


def mass_action_4s(s1: float, s2: float, s3: float, s4: float, k_fwd: float) -> float:
    """Mass-action rate for four substrates."""
    return k_fwd * s1 * s2 * s3 * s4


def reversible_mass_action_keq_1s_1p(
    s1: float, p1: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for one substrate, one product."""
    return kf * (s1 - p1 / keq)


def reversible_mass_action_keq_1s_2p(
    s1: float, p1: float, p2: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for one substrate, two products."""
    return kf * (s1 - p1 * p2 / keq)


def reversible_mass_action_keq_1s_3p(
    s1: float, p1: float, p2: float, p3: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for one substrate, three products."""
    return kf * (s1 - p1 * p2 * p3 / keq)


def reversible_mass_action_keq_2s_1p(
    s1: float, s2: float, p1: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for two substrates, one product."""
    return kf * (s1 * s2 - p1 / keq)


def reversible_mass_action_keq_2s_2p(
    s1: float, s2: float, p1: float, p2: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for two substrates, two products."""
    return kf * (s1 * s2 - p1 * p2 / keq)


def reversible_mass_action_keq_2s_3p(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    p3: float,
    kf: float,
    keq: float,
) -> float:
    """Reversible mass-action rate for two substrates, three products."""
    return kf * (s1 * s2 - p1 * p2 * p3 / keq)


def reversible_mass_action_keq_3s_1p(
    s1: float, s2: float, s3: float, p1: float, kf: float, keq: float
) -> float:
    """Reversible mass-action rate for three substrates, one product."""
    return kf * (s1 * s2 * s3 - p1 / keq)


def reversible_mass_action_keq_3s_2p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    kf: float,
    keq: float,
) -> float:
    """Reversible mass-action rate for three substrates, two products."""
    return kf * (s1 * s2 * s3 - p1 * p2 / keq)


def reversible_mass_action_keq_3s_3p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    kf: float,
    keq: float,
) -> float:
    """Reversible mass-action rate for three substrates, three products."""
    return kf * (s1 * s2 * s3 - p1 * p2 * p3 / keq)


def reversible_mass_action_keq_3s_4p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    kf: float,
    keq: float,
) -> float:
    """Reversible mass-action rate for three substrates, four products."""
    return kf * (s1 * s2 * s3 - p1 * p2 * p3 * p4 / keq)


###############################################################################
# Irreversible Michaelis-Menten type
###############################################################################


def michaelis_menten_1s(
    s: float,
    vmax: float,
    km: float,
) -> float:
    """Irreversible Michaelis-Menten rate for one substrate."""
    return vmax * s / (km + s)


def michaelis_menten_1s_1i(
    s: float,
    i: float,
    vmax: float,
    km: float,
    ki: float,
) -> float:
    """Irreversible Michaelis-Menten rate for one substrate with one inhibitor."""
    return vmax * s / (s + km * (1 + i / ki))


def michaelis_menten_1s_2i(
    s: float,
    i1: float,
    i2: float,
    vmax: float,
    km: float,
    ki1: float,
    ki2: float,
) -> float:
    """Irreversible Michaelis-Menten rate for one substrate with two inhibitors."""
    return vmax * s / (s + km * (1 + i1 / ki1 + i2 / ki2))


def michaelis_menten_1s_3i(
    s: float,
    i1: float,
    i2: float,
    i3: float,
    vmax: float,
    km: float,
    ki1: float,
    ki2: float,
    ki3: float,
) -> float:
    """Irreversible Michaelis-Menten rate for one substrate with three inhibitors."""
    return vmax * s / (s + km * (1 + i1 / ki1 + i2 / ki2 + i3 / ki3))


def michaelis_menten_1s_4i(
    s: float,
    i1: float,
    i2: float,
    i3: float,
    i4: float,
    vmax: float,
    km: float,
    ki1: float,
    ki2: float,
    ki3: float,
    ki4: float,
) -> float:
    """Irreversible Michaelis-Menten rate for one substrate with four inhibitors."""
    return vmax * s / (s + km * (1 + i1 / ki1 + i2 / ki2 + i3 / ki3 + i4 / ki4))


def michaelis_menten_2s(
    s1: float,
    s2: float,
    vmax: float,
    km: float,
) -> float:
    """Irreversible Michaelis-Menten rate for two substrates, shared Km."""
    return vmax * s1 * s2 / (s1 * s2 + km * s1 + km * s2)


def michaelis_menten_2s_1i(
    s1: float,
    s2: float,
    vmax: float,
    km: float,
    i1: float,
    ki1: float,
) -> float:
    """Irreversible Michaelis-Menten rate for two substrates with one inhibitor."""
    km = km * (1 + i1 / ki1)
    return vmax * s1 * s2 / (s1 * s2 + km * s1 + km * s2)


def michaelis_menten_2s_2km(
    s1: float,
    s2: float,
    vmax: float,
    km_s1: float,
    km_s2: float,
) -> float:
    """Irreversible Michaelis-Menten rate for two substrates with individual Km values."""
    return vmax * s1 * s2 / (s1 * s2 + km_s1 * s1 + km_s2 * s2)


def michaelis_menten_3s(
    s1: float,
    s2: float,
    s3: float,
    vmax: float,
    km: float,
) -> float:
    """Irreversible Michaelis-Menten rate for three substrates, shared Km."""
    return (
        vmax
        * s1
        * s2
        * s3
        / (s1 * s2 * s3 + km * s2 * s3 + km * s1 * s3 + km * s1 * s2)
    )


def michaelis_menten_3s_1i(
    s1: float,
    s2: float,
    s3: float,
    vmax: float,
    km: float,
    i1: float,
    ki1: float,
) -> float:
    """Irreversible Michaelis-Menten rate for three substrates with one inhibitor."""
    km = km * (1 + i1 / ki1)
    return (
        vmax
        * s1
        * s2
        * s3
        / (s1 * s2 * s3 + km * s2 * s3 + km * s1 * s3 + km * s1 * s2)
    )


def michaelis_menten_3s_2i(
    s1: float,
    s2: float,
    s3: float,
    vmax: float,
    km: float,
    i1: float,
    ki1: float,
    i2: float,
    ki2: float,
) -> float:
    """Irreversible Michaelis-Menten rate for three substrates with two inhibitors."""
    km = km * (1 + i1 / ki1 + i2 / ki2)
    return (
        vmax
        * s1
        * s2
        * s3
        / (s1 * s2 * s3 + km * s2 * s3 + km * s1 * s3 + km * s1 * s2)
    )


###############################################################################
# Reversible Michaelis-Menten type
###############################################################################


def reversible_michaelis_menten_1s_1p(
    s1: float,
    p1: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for one substrate, one product."""
    return (vmax / kms / (1 + s1 / kms + p1 / kmp)) * (s1 - p1 / keq)


def reversible_michaelis_menten_1s_1p_1i(
    s1: float,
    p1: float,
    i1: float,
    vmax: float,
    kms: float,
    kmp: float,
    ki: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for one substrate, one product, one inhibitor."""
    return (vmax / kms / (1 + s1 / (kms * (1 + i1 / ki)) + p1 / kmp)) * (s1 - p1 / keq)


def reversible_michaelis_menten_1s_2p(
    s1: float,
    p1: float,
    p2: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for one substrate, two products."""
    return (vmax / kms / (1 + s1 / kms + p1 * p2 / kmp)) * (s1 - p1 * p2 / keq)


def reversible_michaelis_menten_1s_3p(
    s1: float,
    p1: float,
    p2: float,
    p3: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for one substrate, three products."""
    return (vmax / kms / (1 + s1 / kms + p1 * p2 * p3 / kmp)) * (
        s1 - p1 * p2 * p3 / keq
    )


def reversible_michaelis_menten_2s_1p(
    s1: float,
    s2: float,
    p1: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for two substrates, one product."""
    return (vmax / kms / (1 + s1 * s2 / kms + p1 / kmp)) * (s1 * s2 - p1 / keq)


def reversible_michaelis_menten_2s_2p(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for two substrates, two products."""
    return (vmax / kms / (1 + s1 * s2 / kms + p1 * p2 / kmp)) * (
        s1 * s2 - p1 * p2 / keq
    )


def reversible_michaelis_menten_2s_3p(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    p3: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for two substrates, three products."""
    return (vmax / kms / (1 + s1 * s2 / kms + p1 * p2 * p3 / kmp)) * (
        s1 * s2 - p1 * p2 * p3 / keq
    )


def reversible_michaelis_menten_2s_4p(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for two substrates, four products."""
    return (vmax / kms / (1 + s1 * s2 / kms + p1 * p2 * p3 * p4 / kmp)) * (
        s1 * s2 - p1 * p2 * p3 * p4 / keq
    )


def reversible_michaelis_menten_3s_1p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, one product."""
    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 / kmp)) * (
        s1 * s2 * s3 - p1 / keq
    )


def reversible_michaelis_menten_3s_2p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, two products."""
    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 * p2 / kmp)) * (
        s1 * s2 * s3 - p1 * p2 / keq
    )


def reversible_michaelis_menten_3s_3p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, three products."""
    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 * p2 * p3 / kmp)) * (
        s1 * s2 * s3 - p1 * p2 * p3 / keq
    )


def reversible_michaelis_menten_3s_3p_1i(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
    i1: float,
    ki1: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, three products, one inhibitor."""
    kms = kms * (1 + i1 / ki1)
    kmp = kmp * (1 + i1 / ki1)

    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 * p2 * p3 / kmp)) * (
        s1 * s2 * s3 - p1 * p2 * p3 / keq
    )


def reversible_michaelis_menten_3s_3p_2i(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
    i1: float,
    ki1: float,
    i2: float,
    ki2: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, three products, two inhibitors."""
    kms = kms * (1 + i1 / ki1 + i2 / ki2)
    kmp = kmp * (1 + i1 / ki1 + i2 / ki2)

    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 * p2 * p3 / kmp)) * (
        s1 * s2 * s3 - p1 * p2 * p3 / keq
    )


def reversible_michaelis_menten_3s_4p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    vmax: float,
    kms: float,
    kmp: float,
    keq: float,
) -> float:
    """Reversible Michaelis-Menten rate for three substrates, four products."""
    return (vmax / kms / (1 + s1 * s2 * s3 / kms + p1 * p2 * p3 * p4 / kmp)) * (
        s1 * s2 * s3 - p1 * p2 * p3 * p4 / keq
    )


def ping_pong_bi_bi(
    s1: float, s2: float, vmax: float, km_s1: float, km_s2: float
) -> float:
    """Ping-pong Bi-Bi kinetics for two-substrate, two-product reactions."""
    return vmax * s1 * s2 / (1 / (km_s1 * km_s2) + s1 / km_s1 + s2 / km_s2 + s1 * s2)


def rapid_equilibrium_1s_1p(
    s1: float,
    p1: float,
    k_re: float,
    q: float,
) -> float:
    """Rapid-equilibrium rate for one substrate, one product."""
    return k_re * (s1 - p1 / q)


def rapid_equilibrium_2s_1p(
    s1: float,
    s2: float,
    p1: float,
    k_re: float,
    q: float,
) -> float:
    """Rapid-equilibrium rate for two substrates, one product."""
    return k_re * (s1 * s2 - p1 / q)


def rapid_equilibrium_2s_2p(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    k_re: float,
    q: float,
) -> float:
    """Rapid-equilibrium rate for two substrates, two products."""
    return k_re * (s1 * s2 - p1 * p2 / q)


def rapid_equilibrium_3s_3p(
    s1: float,
    s2: float,
    s3: float,
    p1: float,
    p2: float,
    p3: float,
    k_re: float,
    q: float,
) -> float:
    """Rapid-equilibrium rate for three substrates, three products."""
    return k_re * (s1 * s2 * s3 - p1 * p2 * p3 / q)


# Misc


def protons_stroma(ph: float) -> float:
    """Convert stromal pH to proton concentration (µmol/L)."""
    return 4000.0 * 10 ** (-ph)

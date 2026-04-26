"""Enzyme building block: b6f."""

import math

from mxlpy import Derived, Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    static,
)


def _four_div_by(x: float) -> float:
    return 4.0 / x


def _keq_cytb6f(
    pH: float,
    F: float,
    E0_PQ: float,
    E0_PC: float,
    pHstroma: float,
    RT: float,
    dG_pH: float,
) -> float:
    DG1 = -2 * F * E0_PQ
    DG2 = -F * E0_PC
    DG = -(DG1 + 2 * dG_pH * pH) + 2 * DG2 + 2 * dG_pH * (pHstroma - pH)
    return math.exp(-DG / RT)


def _b6f(
    PC_ox: float,
    PQ_ox: float,
    PQ_red: float,
    PC_red: float,
    Keq_B6f: float,
    kCytb6f: float,
) -> float:
    return max(
        kCytb6f * (PQ_red * PC_ox**2 - PQ_ox * PC_red**2 / Keq_B6f),
        -kCytb6f,
    )


def _k_b6f(
    pH: float,
    pKreg: float,
    b6f_content: float,
    max_b6f: float,
) -> float:
    pHmod = 1 - (1 / (10 ** (pH - pKreg) + 1))
    b6f_deprot = pHmod * b6f_content
    return b6f_deprot * max_b6f


def _b6f_2024(
    PC: float,
    PCred: float,
    PQ: float,
    PQred: float,
    k_b6f: float,
    Keq_cytb6f: float,
) -> float:
    k_b6f_reverse = k_b6f / Keq_cytb6f
    f_PQH2 = PQred / (
        PQred + PQ
    )  # want to keep the rates in terms of fraction of PQHs, not total number
    f_PQ = 1 - f_PQH2
    return f_PQH2 * PC * k_b6f - f_PQ * PCred * k_b6f_reverse


def add_b6f(
    model: Model,
    *,
    rxn: str | None = None,
    ph_stroma: str | None = None,
    ph_lumen: str | None = None,
    h_lumen: str | None = None,
    pc_ox: str | None = None,
    pq_ox: str | None = None,
    pq_red: str | None = None,
    pc_red: str | None = None,
    bh: str | None = None,
) -> Model:
    """Add b6f to model."""
    rxn = default_name(rxn, n.b6f)
    ph_stroma = default_name(ph_stroma, n.ph)
    ph_lumen = default_name(ph_lumen, lambda: n.ph("_lumen"))
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    pc_ox = default_name(pc_ox, n.pc_ox)
    pq_ox = default_name(pq_ox, n.pq_ox)
    pq_red = default_name(pq_red, n.pq_red)
    pc_red = default_name(pc_red, n.pc_red)

    bh = static(model, "bH", 100.0) if bh is None else bh

    model.add_parameter(n.kcat(rxn), 2.5)
    model.add_derived(
        name=n.keq(rxn),
        fn=_keq_cytb6f,
        args=[
            ph_lumen,
            "F",
            "E^0_PQ",
            "E^0_PC",
            ph_stroma,
            "RT",
            "dG_pH",
        ],
    )
    model.add_reaction(
        name=rxn,
        fn=_b6f,
        stoichiometry={
            pc_ox: -2,
            pq_ox: 1,
            h_lumen: Derived(fn=_four_div_by, args=[bh]),
        },
        args=[
            pc_ox,
            pq_ox,
            pq_red,
            pc_red,
            n.keq(rxn),
            n.kcat(rxn),
        ],
    )
    return model


def add_b6f_2024(
    model: Model,
    *,
    rxn: str | None = None,
    ph_stroma: str | None = None,
    ph_lumen: str | None = None,
    h_lumen: str | None = None,
    pc_ox: str | None = None,
    pq_ox: str | None = None,
    pq_red: str | None = None,
    pc_red: str | None = None,
) -> Model:
    """Add b6f 2024 to model."""
    rxn = default_name(rxn, n.b6f)
    ph_stroma = default_name(ph_stroma, n.ph)
    ph_lumen = default_name(ph_lumen, lambda: n.ph("_lumen"))
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))
    pc_ox = default_name(pc_ox, n.pc_ox)
    pq_ox = default_name(pq_ox, n.pq_ox)
    pq_red = default_name(pq_red, n.pq_red)
    pc_red = default_name(pc_red, n.pc_red)

    model.add_parameter(b6f_content := "b6f_content", 1)
    model.add_parameter(max_b6f := "max_b6f", 500)
    model.add_parameter(pKreg := "pKreg", 6.5)

    model.add_derived(
        name=n.keq(rxn),
        fn=_keq_cytb6f,
        args=[
            ph_lumen,
            "F",
            "E^0_PQ",
            "E^0_PC",
            ph_stroma,
            "RT",
            "dG_pH",
        ],
    )

    model.add_derived(
        name=n.keq(rxn + "_dyn"),
        fn=_k_b6f,
        args=[
            ph_lumen,
            pKreg,
            b6f_content,
            max_b6f,
        ],
    )

    model.add_reaction(
        name=rxn,
        fn=_b6f_2024,
        stoichiometry={
            pc_ox: -2,
            pq_ox: 1,
            h_lumen: Derived(fn=_four_div_by, args=["bH"]),
        },
        args=[
            pc_ox,
            pc_red,
            pq_ox,
            pq_red,
            n.keq(rxn + "_dyn"),
            n.keq(rxn),
        ],
    )

    return model

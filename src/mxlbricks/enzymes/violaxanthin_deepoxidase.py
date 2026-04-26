"""Violaxanthin Deepoxidase (lumen).

Violaxanthin + Ascorbate -> Antheraxanthin + Dehydroascorbate + H2O
Antheraxanthin + Ascorbate -> Zeaxanthin + Dehydroascorbate + H2O.
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import protons_stroma
from mxlbricks.utils import (
    default_name,
    static,
)


def _rate_protonation_hill(
    Vx: float,
    H: float,
    k_fwd: float,
    nH: float,
    kphSat: float,
) -> float:
    return k_fwd * (H**nH / (H**nH + protons_stroma(kphSat) ** nH)) * Vx  # type: ignore


def add_violaxanthin_epoxidase(
    model: Model,
    *,
    rxn: str | None = None,
    vx: str | None = None,
    h_lumen: str | None = None,
    kf: str | None = None,
    kh_zx: str | None = None,
    kphsat: str | None = None,
) -> Model:
    """Add violaxanthin epoxidase (Violaxanthin Deepoxidase (lumen)) to model."""
    rxn = default_name(rxn, n.violaxanthin_deepoxidase)
    vx = default_name(vx, n.vx)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))

    model.add_reaction(
        name=rxn,
        fn=_rate_protonation_hill,
        stoichiometry={
            vx: -1,
        },
        args=[
            vx,
            h_lumen,
            static(model, n.kf(rxn), 0.0024) if kf is None else kf,
            static(model, n.kh(rxn), 5.0) if kh_zx is None else kh_zx,
            static(model, n.ksat(rxn), 5.8) if kphsat is None else kphsat,
        ],
    )
    return model

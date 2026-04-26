"""name.

EC FIXME

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import protons_stroma
from mxlbricks.utils import (
    default_name,
    filter_stoichiometry,
    static,
)


def _protonation_hill(
    vx: float,
    h: float,
    nh: float,
    k_fwd: float,
    k_ph_sat: float,
) -> float:
    return k_fwd * (h**nh / (h**nh + protons_stroma(k_ph_sat) ** nh)) * vx  # type: ignore


def add_lhc_protonation(
    model: Model,
    *,
    rxn: str | None = None,
    psbs_de: str | None = None,
    psbs_pr: str | None = None,
    h_lumen: str | None = None,
    kf: str | None = None,
    kh_lhc: str | None = None,
    k_ph_sat: str | None = None,
) -> Model:
    """Add lhc protonation (name) to model."""
    rxn = default_name(rxn, n.lhc_protonation)
    psbs_de = default_name(psbs_de, n.psbs_de)
    psbs_pr = default_name(psbs_pr, n.psbs_pr)
    h_lumen = default_name(h_lumen, lambda: n.h("_lumen"))

    model.add_reaction(
        name=rxn,
        fn=_protonation_hill,
        stoichiometry=filter_stoichiometry(
            model,
            {
                psbs_de: -1,
                psbs_pr: 1,
            },
        ),
        args=[
            psbs_de,
            h_lumen,
            static(model, n.kh(rxn), 3.0) if kh_lhc is None else kh_lhc,
            static(model, n.kf(rxn), 0.0096) if kf is None else kf,
            static(model, n.ksat(rxn), 5.8) if k_ph_sat is None else k_ph_sat,
        ],
    )
    return model

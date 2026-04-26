"""name.

EC 1.8.1.7

glutathione + NADP <=> glutathion-disulfide + NADPH + H+

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_km,
    default_name,
    default_vmax,
)


def _rate_glutathion_reductase(
    nadph: float,
    gssg: float,
    vmax: float,
    km_nadph: float,
    km_gssg: float,
) -> float:
    nom = vmax * nadph * gssg
    denom = km_nadph * gssg + km_gssg * nadph + nadph * gssg + km_nadph * km_gssg
    return nom / denom


def add_glutathion_reductase_irrev(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    glutathion_ox: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_gssg: str | None = None,
    km_nadph: str | None = None,
) -> Model:
    """Add glutathion reductase irrev (name) to model."""
    rxn = default_name(rxn, n.glutathion_reductase)
    nadph = default_name(nadph, n.nadph)
    glutathion_ox = default_name(glutathion_ox, n.glutathion_ox)

    model.add_reaction(
        name=rxn,
        fn=_rate_glutathion_reductase,
        stoichiometry={
            nadph: -1,
            glutathion_ox: -1,
        },
        args=[
            nadph,
            glutathion_ox,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.4e-3,  # Source
                kcat_value=595,  # Source
            ),
            default_km(model, par=km_nadph, rxn=rxn, subs=nadph, value=3e-3),
            default_km(model, par=km_gssg, rxn=rxn, subs=glutathion_ox, value=2e-1),
        ],
    )
    return model

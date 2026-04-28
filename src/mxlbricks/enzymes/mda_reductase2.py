"""EC 1.6.5.4.

NADH + Proton + 2 Monodehydroascorbate <=> NAD + 2 ascorbate.


Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_km,
    default_name,
    default_vmax,
)


def _rate_mda_reductase2(
    nadph: float,
    mda: float,
    vmax: float,
    km_nadph: float,
    km_mda: float,
) -> float:
    """Compare Valero et al. 2016."""
    nom = vmax * nadph * mda
    denom = km_nadph * mda + km_mda * nadph + nadph * mda + km_nadph * km_mda
    return nom / denom


def add_mda_reductase2(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    mda: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_mda: str | None = None,
    km_nadph: str | None = None,
) -> Model:
    """Add mda reductase2 (EC 1) to model."""
    rxn = default_name(rxn, n.mda_reductase2)
    nadph = default_name(nadph, n.nadph)
    mda = default_name(mda, n.mda)

    model.add_reaction(
        name=rxn,
        fn=_rate_mda_reductase2,
        stoichiometry={
            nadph: -1,
            mda: -2,
        },
        args=[
            nadph,
            mda,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=2e-3,  # Source
                kcat_value=1080000 / (60 * 60),  # Source
            ),
            default_km(model, par=km_nadph, rxn=rxn, subs=nadph, value=23e-3),
            default_km(model, par=km_mda, rxn=rxn, subs=mda, value=1.4e-3),
        ],
    )
    return model

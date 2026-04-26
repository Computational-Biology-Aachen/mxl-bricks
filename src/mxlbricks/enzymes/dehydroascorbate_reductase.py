"""dehydroascorbate_reductase, DHAR.

EC FIXME

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    default_vmax,
    static,
)


def _rate_dhar(
    dha: float,
    gsh: float,
    vmax: float,
    km_dha: float,
    km_gsh: float,
    k: float,
) -> float:
    nom = vmax * dha * gsh
    denom = k + km_dha * gsh + km_gsh * dha + dha * gsh
    return nom / denom


def add_dehydroascorbate_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    dha: str | None = None,
    glutathion_ox: str | None = None,
    glutathion_red: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_dha: str | None = None,
    km_gsh: str | None = None,
) -> Model:
    """Add dehydroascorbate reductase (dehydroascorbate_reductase, DHAR) to model."""
    rxn = default_name(rxn, n.dehydroascorbate_reductase)
    dha = default_name(dha, n.dha)
    glutathion_ox = default_name(glutathion_ox, n.glutathion_ox)
    glutathion_red = default_name(glutathion_red, n.glutathion_red)

    km_dha = (
        static(model, n.km(rxn, n.dha()), 70e-3) if km_dha is None else km_dha
    )  # FIXME: source
    km_gsh = (
        static(model, n.km(rxn, n.glutathion_red()), 2.5e3 * 1e-3)
        if km_gsh is None
        else km_gsh
    )  # FIXME: source

    model.add_parameter("K", 5e5 * (1e-3) ** 2)

    model.add_reaction(
        name=rxn,
        fn=_rate_dhar,
        stoichiometry={
            dha: -1,
            glutathion_ox: 1,
        },
        args=[
            dha,
            glutathion_red,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.7e-3,  # Source
                kcat_value=142,  # Source
            ),
            km_dha,
            km_gsh,
            "K",
        ],
    )
    return model

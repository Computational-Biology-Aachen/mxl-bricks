"""ascorbate peroxidase.

EC FIXME

Equilibrator
Ascorbate(aq) + H2O2(aq) ⇌ Dehydroascorbate(aq) + 2 H2O(l)

Ascorbate(aq) + H2O2(aq) + 0.5 NADH(aq) ⇌ Monodehydroascorbate(aq) + 2 H2O(l) + 0.5 NAD (aq)


"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    filter_stoichiometry,
)


def _rate_ascorbate_peroxidase(
    A: float,
    H: float,
    kf1: float,
    kr1: float,
    kf2: float,
    kr2: float,
    kf3: float,
    kf4: float,
    kr4: float,
    kf5: float,
    XT: float,
) -> float:
    """Lumped reaction of ascorbate peroxidase.

    the cycle stretched to a linear chain with
    two steps producing the MDA
    two steps releasing ASC
    and one step producing hydrogen peroxide.
    """
    nom = A * H * XT
    denom = (
        A * H * (1 / kf3 + 1 / kf5)
        + A / kf1
        + H / kf4
        + H * kr4 / (kf4 * kf5)
        + H / kf2
        + H * kr2 / (kf2 * kf3)
        + kr1 / (kf1 * kf2)
        + kr1 * kr2 / (kf1 * kf2 * kf3)
    )
    return nom / denom


def add_ascorbate_peroxidase(
    model: Model,
    *,
    rxn: str | None = None,
    s1: str | None = None,
    s2: str | None = None,
    p1: str | None = None,
) -> Model:
    """Add ascorbate peroxidase (ascorbate peroxidase) to model."""
    rxn = default_name(rxn, n.ascorbate_peroxidase)
    s1 = default_name(s1, n.ascorbate)
    s2 = default_name(s2, n.h2o2)
    p1 = default_name(p1, n.mda)

    model.add_parameters(
        {
            "kf1": 10000.0,
            "kr1": 220.0,
            "kf2": 10000.0,
            "kr2": 4000.0,
            "kf3": 2510.0,
            "kf4": 10000.0,
            "kr4": 4000.0,
            "kf5": 2510.0,
            "XT": 0.07,  # according to Valero
        }
    )
    model.add_reaction(
        name=rxn,
        fn=_rate_ascorbate_peroxidase,
        stoichiometry=filter_stoichiometry(
            model,
            {
                s1: -1,
                s2: -1,
                p1: 2,
            },
        ),
        args=[
            s1,
            s2,
            "kf1",
            "kr1",
            "kf2",
            "kr2",
            "kf3",
            "kf4",
            "kr4",
            "kf5",
            "XT",
        ],
    )
    return model

"""GAP + F6P <=> E4P + X5P.

EC 2.2.1.1

Equilibrator
D-Glyceraldehyde 3-phosphate(aq) + D-Fructose 6-phosphate(aq)
    ⇌ D-Xylulose 5-phosphate(aq) + D-Erythrose 4-phosphate(aq)
Keq = 0.02 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_2s_2p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_transketolase_x5p_e4p_f6p_gap(
    model: Model,
    *,
    rxn: str | None = None,
    gap: str | None = None,
    f6p: str | None = None,
    e4p: str | None = None,
    x5p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add transketolase x5p e4p f6p gap (GAP + F6P <=> E4P + X5P) to model."""
    rxn = default_name(rxn, n.transketolase_gap_f6p)
    gap = default_name(gap, n.gap)
    f6p = default_name(f6p, n.f6p)
    e4p = default_name(e4p, n.e4p)
    x5p = default_name(x5p, n.x5p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_2p,
        stoichiometry={
            gap: -1,
            f6p: -1,
            e4p: 1,
            x5p: 1,
        },
        args=[
            gap,
            f6p,
            e4p,
            x5p,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=0.084),
        ],
    )

    return model

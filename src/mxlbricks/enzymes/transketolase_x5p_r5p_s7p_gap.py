"""GAP + S7P <=> R5P + X5P.

EC 2.2.1.1

Equilibrator
D-Glyceraldehyde 3-phosphate(aq) + Sedoheptulose 7-phosphate(aq)
    ⇌ D-Ribose 5-phosphate(aq) + D-Xylulose 5-phosphate(aq)
Keq = 0.2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_2s_2p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
)


def add_transketolase_x5p_r5p_s7p_gap(
    model: Model,
    *,
    rxn: str | None = None,
    gap: str | None = None,
    s7p: str | None = None,
    r5p: str | None = None,
    x5p: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add transketolase x5p r5p s7p gap (GAP + S7P <=> R5P + X5P) to model."""
    rxn = default_name(rxn, n.transketolase_gap_s7p)
    gap = default_name(gap, n.gap)
    s7p = default_name(s7p, n.s7p)
    r5p = default_name(r5p, n.r5p)
    x5p = default_name(x5p, n.x5p)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_2s_2p,
        stoichiometry={
            gap: -1,
            s7p: -1,
            r5p: 1,
            x5p: 1,
        },
        args=[
            gap,
            s7p,
            r5p,
            x5p,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=0.85),
        ],
    )
    return model

"""Glyceraldehyde 3-phosphate dehydrogenase (GADPH).

EC 1.2.1.13

Equilibrator
    NADPH(aq) + 3-Phospho-D-glyceroyl phosphate(aq)
    ⇌ NADP (aq) + Orthophosphate(aq) + D-Glyceraldehyde 3-phosphate(aq)
    Keq = 2 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)

FIXME: Poolman uses H+ in the description. Why?
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import rapid_equilibrium_3s_3p
from mxlbricks.utils import (
    default_keq,
    default_kre,
    default_name,
    filter_stoichiometry,
)


def add_gadph(
    model: Model,
    *,
    rxn: str | None = None,
    bpga: str | None = None,
    nadph: str | None = None,
    h: str | None = None,
    gap: str | None = None,
    nadp: str | None = None,
    pi: str | None = None,
    kre: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add gadph (Glyceraldehyde 3-phosphate dehydrogenase (GADPH)) to model."""
    rxn = default_name(rxn, n.gadph)
    bpga = default_name(bpga, n.bpga)
    nadph = default_name(nadph, n.nadph)
    h = default_name(h, n.h)
    gap = default_name(gap, n.gap)
    nadp = default_name(nadp, n.nadp)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=rapid_equilibrium_3s_3p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadph: -1.0,
                bpga: -1.0,
                nadp: 1.0,
                pi: 1.0,
                gap: 1.0,
            },
        ),
        args=[
            bpga,
            nadph,
            h,
            gap,
            nadp,
            pi,
            default_kre(model, par=kre, rxn=rxn, value=800000000.0),
            default_keq(model, rxn=rxn, par=keq, value=16000000.0),
        ],
    )
    return model

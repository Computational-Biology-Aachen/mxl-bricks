"""EC 5.4.2.1

Equilibrator
2-Phospho-D-glycerate(aq) ⇌ 3-Phospho-D-glycerate(aq)
Keq = 6 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, reversible_michaelis_menten_1s_1p
from mxlbricks.utils import filter_stoichiometry, static

if TYPE_CHECKING:
    from mxlpy import Model

ENZYME = n.phosphoglycerate_mutase()


def add_phosphoglycerate_mutase(
    model: Model,
    compartment: str = "",
    *,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    kms = static(model, n.kms(ENZYME), 0.1) if kms is None else kms  # FIXME: source
    kmp = static(model, n.kmp(ENZYME), 1.0) if kmp is None else kmp  # FIXME: source
    kcat = static(model, n.kcat(ENZYME), 1.0) if kcat is None else kcat  # FIXME: source
    e0 = static(model, n.e0(ENZYME), 1.0) if e0 is None else e0  # FIXME: source
    keq = static(model, n.keq(ENZYME), 6.0) if keq is None else keq  # FIXME: source
    model.add_derived(vmax := n.vmax(ENZYME), fn=mass_action_1s, args=[kcat, e0])

    model.add_reaction(
        name=ENZYME,
        fn=reversible_michaelis_menten_1s_1p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                n.pga2(compartment): -1.0,
                n.pga(compartment): 1.0,
            },
        ),
        args=[
            n.pga2(compartment),
            n.pga(compartment),
            vmax,
            kms,
            kmp,
            keq,
        ],
    )
    return model

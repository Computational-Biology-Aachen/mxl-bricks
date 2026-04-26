"""EC 3.1.2.10.

Metacyc:
FORMYL_COA_m + WATER_m <=> CO-A_m + FORMATE_m + PROTON_m

Equilibrator
Formyl-CoA(aq) + Water(l) ⇌ CoA(aq) + Formate(aq)
Keq = 4e1 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, reversible_michaelis_menten_1s_2p
from mxlbricks.utils import (
    default_name,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_thioesterase(
    model: Model,
    *,
    rxn: str | None = None,
    formyl_coa: str | None = None,
    coa: str | None = None,
    formate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add thioesterase (EC 3) to model."""
    rxn = default_name(rxn, n.thioesterase)
    formyl_coa = default_name(formyl_coa, n.formyl_coa)
    coa = default_name(coa, n.coa)
    formate = default_name(formate, n.formate)

    kms = static(model, n.kms(rxn), 0.045) if kms is None else kms  # FIXME: source
    kmp = static(model, n.kmp(rxn), 1.0) if kmp is None else kmp  # FIXME: source
    kcat = static(model, n.kcat(rxn), 1.0) if kcat is None else kcat  # FIXME: source
    e0 = static(model, n.e0(rxn), 1.0) if e0 is None else e0  # FIXME: source
    keq = static(model, n.keq(rxn), 40.0) if keq is None else keq  # FIXME: source
    model.add_derived(vmax := n.vmax(rxn), fn=mass_action_1s, args=[kcat, e0])

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                formyl_coa: -1,
                coa: 1,
                formate: 1,
            },
        ),
        args=[
            formyl_coa,
            coa,
            formate,
            vmax,
            kms,
            kmp,
            keq,
        ],
    )

    return model

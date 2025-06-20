"""glyoxylate carboligase == tartronate-semialdehyde synthase

EC 4.1.1.47

Equilibrator
2 Glyoxylate + H2O <=> Tartronate semialdehyde + CO2(total)
Keq = 1.6e4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s, reversible_michaelis_menten_1s_2p
from mxlbricks.utils import filter_stoichiometry, static

if TYPE_CHECKING:
    from mxlpy import Model

ENZYME = n.glyoxylate_carboligase()


def add_glyoxylate_carboligase(
    model: Model,
    compartment: str = "",
    *,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    kms = static(model, n.kms(ENZYME), 0.9) if kms is None else kms  # FIXME: source
    kmp = static(model, n.kmp(ENZYME), 1.0) if kmp is None else kmp  # FIXME: source
    kcat = (
        static(model, n.kcat(ENZYME), 18.9) if kcat is None else kcat
    )  # FIXME: source
    e0 = static(model, n.e0(ENZYME), 1.0) if e0 is None else e0  # FIXME: source
    keq = static(model, n.keq(ENZYME), 1.6e4) if keq is None else keq  # FIXME: source
    model.add_derived(vmax := n.vmax(ENZYME), fn=mass_action_1s, args=[kcat, e0])

    model.add_reaction(
        name=ENZYME,
        fn=reversible_michaelis_menten_1s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                n.glyoxylate(compartment): -2,
                n.tartronate_semialdehyde(compartment): 1,
                n.co2(compartment): 1,
            },
        ),
        args=[
            n.glyoxylate(compartment),
            n.tartronate_semialdehyde(compartment),
            n.co2(compartment),
            vmax,
            kms,
            kmp,
            keq,
        ],
    )
    return model

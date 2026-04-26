"""Enzyme building block: oxidative phosphorylation."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import reversible_mass_action_keq_2s_2p
from mxlbricks.utils import (
    default_keq,
    default_kf,
    default_name,
    filter_stoichiometry,
)


def add_oxidative_phosphorylation(
    model: Model,
    *,
    rxn: str | None = None,
    nadph: str | None = None,
    adp: str | None = None,
    nadp: str | None = None,
    atp: str | None = None,
    kf: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add oxidative phosphorylation to model."""
    rxn = default_name(rxn, n.oxidative_phosphorylation)
    nadph = default_name(nadph, n.nadph)
    adp = default_name(adp, n.adp)
    nadp = default_name(nadp, n.nadp)
    atp = default_name(atp, n.atp)

    model.add_reaction(
        rxn,
        reversible_mass_action_keq_2s_2p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                nadph: -1,
                adp: -1,
                nadp: 1,
                atp: 1,
            },
        ),
        args=[
            nadph,
            adp,
            nadp,
            atp,
            default_kf(model, rxn=rxn, par=kf, value=1),
            default_keq(model, rxn=rxn, par=keq, value=3 / 2),
        ],
    )
    return model

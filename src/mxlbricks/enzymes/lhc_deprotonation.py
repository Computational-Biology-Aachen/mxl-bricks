"""Enzyme building block: lhc deprotonation."""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_kf,
    default_name,
    filter_stoichiometry,
)


def add_lhc_deprotonation(
    model: Model,
    *,
    rxn: str | None = None,
    psbs_pr: str | None = None,
    psbs_de: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add lhc deprotonation to model."""
    rxn = default_name(rxn, n.lhc_deprotonation)
    psbs_pr = default_name(psbs_pr, n.psbs_pr)
    psbs_de = default_name(psbs_de, n.psbs_de)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                psbs_pr: -1,
                psbs_de: 1,
            },
        ),
        args=[
            psbs_pr,
            default_kf(model, rxn=rxn, par=kf, value=0.0096),
        ],
    )
    return model

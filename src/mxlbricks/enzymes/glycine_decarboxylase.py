"""glycine decarboxylase.

2 Glycine + NAD + 2 H2O ⇌ Serine + NH3 + NADH + CO2

Equilibrator
2 Glycine(aq) + NAD(aq) + 2 H2O(l) ⇌ Serine(aq) + NH3(aq) + NADH(aq) + CO2(total)
Keq = 2.4e-4 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    michaelis_menten_1s,
    michaelis_menten_2s,
    reversible_michaelis_menten_2s_4p,
)
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_glycine_decarboxylase_yokota(
    model: Model,
    *,
    rxn: str | None = None,
    glycine: str | None = None,
    serine: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycine decarboxylase yokota (glycine decarboxylase) to model."""
    rxn = default_name(rxn, n.glycine_decarboxylase)

    glycine = default_name(glycine, n.glycine)
    serine = default_name(serine, n.serine)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycine: -2.0,
                serine: 1.0,
            },
        ),
        args=[
            glycine,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=100.0,  # Source
            ),
            default_kms(model, par=kms, rxn=rxn, value=6.0),
        ],
    )
    return model


def add_glycine_decarboxylase_irreversible(
    model: Model,
    *,
    rxn: str | None = None,
    glycine: str | None = None,
    nad: str | None = None,
    serine: str | None = None,
    nh4: str | None = None,
    nadh: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
) -> Model:
    """Add glycine decarboxylase irreversible (glycine decarboxylase) to model."""
    rxn = default_name(rxn, n.glycine_decarboxylase)
    glycine = default_name(glycine, n.glycine)
    nad = default_name(nad, n.nad)
    serine = default_name(serine, n.serine)
    nh4 = default_name(nh4, n.nh4)
    nadh = default_name(nadh, n.nadh)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=michaelis_menten_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycine: -2.0,
                nad: -1.0,
                serine: 1.0,
                nh4: 1.0,
                nadh: 1.0,
                co2: 1.0,
            },
        ),
        args=[
            glycine,
            nad,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=100.0,  # Source
            ),
            static(model, n.kms(rxn), 6.0) if kms is None else kms,
        ],
    )

    return model


def add_glycine_decarboxylase(
    model: Model,
    *,
    rxn: str | None = None,
    glycine: str | None = None,
    nad: str | None = None,
    serine: str | None = None,
    nh4: str | None = None,
    nadh: str | None = None,
    co2: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add glycine decarboxylase (glycine decarboxylase) to model."""
    rxn = default_name(rxn, n.glycine_decarboxylase)
    glycine = default_name(glycine, n.glycine)
    nad = default_name(nad, n.nad)
    serine = default_name(serine, n.serine)
    nh4 = default_name(nh4, n.nh4)
    nadh = default_name(nadh, n.nadh)
    co2 = default_name(co2, n.co2)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_4p,
        stoichiometry=filter_stoichiometry(
            model,
            {
                glycine: -2.0,
                nad: -1.0,
                serine: 1.0,
                nh4: 1.0,
                nadh: 1.0,
                co2: 1.0,
            },
        ),
        args=[
            glycine,
            nad,
            serine,
            nh4,
            nadh,
            co2,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=100.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=6.0),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=0.00024),
        ],
    )
    return model

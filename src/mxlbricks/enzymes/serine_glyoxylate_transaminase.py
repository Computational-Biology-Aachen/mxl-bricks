"""serine glyoxylate transaminase.

Glyoxylate + L-Serine <=> Glycine + Hydroxypyruvate

EC 2.6.1.45

Equilibrator
Glyoxylate(aq) + Serine(aq) ⇌ Glycine(aq) + Hydroxypyruvate(aq)
Keq = 6 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import (
    ping_pong_bi_bi,
    reversible_michaelis_menten_2s_2p,
)
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_serine_glyoxylate_transaminase_irreversible(
    model: Model,
    *,
    rxn: str | None = None,
    glyoxylate: str | None = None,
    serine: str | None = None,
    glycine: str | None = None,
    hydroxypyruvate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_gox: str | None = None,
    km_ser: str | None = None,
) -> Model:
    """Add serine glyoxylate transaminase irreversible (serine glyoxylate transaminase) to model."""
    rxn = default_name(rxn, n.serine_glyoxylate_transaminase)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    serine = default_name(serine, n.serine)
    glycine = default_name(glycine, n.glycine)
    hydroxypyruvate = default_name(hydroxypyruvate, n.hydroxypyruvate)

    model.add_reaction(
        name=rxn,
        fn=ping_pong_bi_bi,
        stoichiometry={
            glyoxylate: -1.0,
            serine: -1.0,
            glycine: 1.0,
            hydroxypyruvate: 1.0,
        },
        args=[
            glyoxylate,
            serine,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=159.0,  # Source
            ),
            static(model, n.km(rxn, n.glyoxylate()), 0.15)
            if km_gox is None
            else km_gox,
            static(model, n.km(rxn, n.serine()), 2.72) if km_ser is None else km_ser,
        ],
    )

    return model


def add_serine_glyoxylate_transaminase(
    model: Model,
    *,
    rxn: str | None = None,
    glyoxylate: str | None = None,
    serine: str | None = None,
    glycine: str | None = None,
    hydroxypyruvate: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_gox: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add serine glyoxylate transaminase (serine glyoxylate transaminase) to model."""
    rxn = default_name(rxn, n.serine_glyoxylate_transaminase)
    glyoxylate = default_name(glyoxylate, n.glyoxylate)
    serine = default_name(serine, n.serine)
    glycine = default_name(glycine, n.glycine)
    hydroxypyruvate = default_name(hydroxypyruvate, n.hydroxypyruvate)

    # FIXME: kms2 missing
    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_2s_2p,
        stoichiometry={
            glyoxylate: -1.0,
            serine: -1.0,
            glycine: 1.0,
            hydroxypyruvate: 1.0,
        },
        args=[
            glyoxylate,
            serine,
            glycine,
            hydroxypyruvate,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=159.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=km_gox, value=0.15),
            default_kmp(model, rxn=rxn, par=kmp, value=2.72),
            default_keq(model, rxn=rxn, par=keq, value=6.0),
        ],
    )
    return model

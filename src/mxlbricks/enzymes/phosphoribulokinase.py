"""phosphoribulokinase.

EC 2.7.1.19

Equilibrator
    ATP(aq) + D-Ribulose 5-phosphate(aq) ⇌ ADP(aq) + D-Ribulose 1,5-bisphosphate(aq)
    Keq = 1e5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.utils import (
    default_name,
    default_vmax,
    filter_stoichiometry,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def _rate_prk(
    ru5p: float,
    atp: float,
    pi: float,
    pga: float,
    rubp: float,
    adp: float,
    v13: float,
    km131: float,
    km132: float,
    ki131: float,
    ki132: float,
    ki133: float,
    ki134: float,
    ki135: float,
) -> float:
    return (
        v13
        * ru5p
        * atp
        / (
            (ru5p + km131 * (1 + pga / ki131 + rubp / ki132 + pi / ki133))
            * (atp * (1 + adp / ki134) + km132 * (1 + adp / ki135))
        )
    )


def add_phosphoribulokinase(
    model: Model,
    *,
    rxn: str | None = None,
    ru5p: str | None = None,
    atp: str | None = None,
    pi: str | None = None,
    pga: str | None = None,
    rubp: str | None = None,
    adp: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    km_ru5p: str | None = None,
    km_atp: str | None = None,
    ki1: str | None = None,
    ki2: str | None = None,
    ki3: str | None = None,
    ki4: str | None = None,
    ki5: str | None = None,
) -> Model:
    """Add phosphoribulokinase (phosphoribulokinase) to model."""
    rxn = default_name(rxn, n.phosphoribulokinase)
    ru5p = default_name(ru5p, n.ru5p)
    atp = default_name(atp, n.atp)
    pi = default_name(pi, n.pi)
    pga = default_name(pga, n.pga)
    rubp = default_name(rubp, n.rubp)
    adp = default_name(adp, n.adp)

    model.add_reaction(
        name=rxn,
        fn=_rate_prk,
        stoichiometry=filter_stoichiometry(
            model,
            {
                ru5p: -1.0,
                atp: -1.0,
                rubp: 1.0,
                adp: 1.0,
            },
        ),
        args=[
            ru5p,
            atp,
            pi,
            pga,
            rubp,
            adp,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=0.9999 * 8,  # Source
            ),
            static(model, n.km(rxn, n.ru5p()), 0.05) if km_ru5p is None else km_ru5p,
            static(model, n.km(rxn, n.atp()), 0.05) if km_atp is None else km_atp,
            static(model, n.ki(rxn, n.pga()), 2.0) if ki1 is None else ki1,
            static(model, n.ki(rxn, n.rubp()), 0.7) if ki2 is None else ki2,
            static(model, n.ki(rxn, n.pi()), 4.0) if ki3 is None else ki3,
            static(model, n.ki(rxn, "4"), 2.5) if ki4 is None else ki4,
            static(model, n.ki(rxn, "5"), 0.4) if ki5 is None else ki5,
        ],
    )
    return model

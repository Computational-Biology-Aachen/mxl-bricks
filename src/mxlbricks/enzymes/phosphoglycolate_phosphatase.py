"""phosphoglycolate phosphatase, EC 3.1.3.18.

H2O(chl) + PGO(chl) <=> Orthophosphate(chl) + Glycolate(chl)

Equilibrator
H2O(l) + PGO(aq) ⇌ Orthophosphate(aq) + Glycolate(aq)
Keq = 3.1e5 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_1s_1p_1i, value
from mxlbricks.utils import (
    default_keq,
    default_kf,
    default_kis,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
    filter_stoichiometry,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_phosphoglycolate_influx(
    model: Model,
    *,
    rxn: str | None = None,
    glycolate: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add phosphoglycolate influx (phosphoglycolate phosphatase, EC 3) to model."""
    rxn = default_name(rxn, n.phosphoglycolate_phosphatase)
    glycolate = default_name(glycolate, n.glycolate)

    model.add_reaction(
        name=rxn,
        fn=value,
        stoichiometry={
            glycolate: 1,
        },
        args=[
            default_kf(model, par=kf, rxn=rxn, value=60.0),
        ],
    )
    return model


def add_phosphoglycolate_phosphatase(
    model: Model,
    *,
    rxn: str | None = None,
    pgo: str | None = None,
    glycolate: str | None = None,
    pi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
    ki_pi: str | None = None,
) -> Model:
    """Add phosphoglycolate phosphatase (phosphoglycolate phosphatase, EC 3) to model."""
    rxn = default_name(rxn, n.phosphoglycolate_phosphatase)
    pgo = default_name(pgo, n.pgo)
    glycolate = default_name(glycolate, n.glycolate)
    pi = default_name(pi, n.pi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_1s_1p_1i,
        stoichiometry=filter_stoichiometry(
            model,
            {
                pgo: -1.0,
                glycolate: 1.0,
                pi: 1.0,
            },
        ),
        args=[
            pgo,
            glycolate,
            pi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=292.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.029),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_kis(model, rxn=rxn, par=ki_pi, substrate=pi, value=12.0),
            default_keq(model, rxn=rxn, par=keq, value=310000.0),
        ],
    )
    return model

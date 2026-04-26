"""name.

EC 2.7.9.1

Equilibrator
Pyruvate(aq) + ATP(aq) + Orthophosphate(aq) ⇌ PEP(aq) + AMP(aq) + Diphosphate(aq)
Keq = 9.6e-3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import reversible_michaelis_menten_3s_3p
from mxlbricks.utils import (
    default_keq,
    default_kmp,
    default_kms,
    default_name,
    default_vmax,
)

if TYPE_CHECKING:
    from mxlpy import Model


def add_pyruvate_phosphate_dikinase(
    model: Model,
    *,
    rxn: str | None = None,
    pyruvate: str | None = None,
    atp: str | None = None,
    pi: str | None = None,
    pep: str | None = None,
    amp: str | None = None,
    ppi: str | None = None,
    kcat: str | None = None,
    e0: str | None = None,
    kms: str | None = None,
    kmp: str | None = None,
    keq: str | None = None,
) -> Model:
    """Add pyruvate phosphate dikinase (name) to model."""
    rxn = default_name(rxn, n.pyruvate_phosphate_dikinase)
    pyruvate = default_name(pyruvate, n.pyruvate)
    atp = default_name(atp, n.atp)
    pi = default_name(pi, n.pi)
    pep = default_name(pep, n.pep)
    amp = default_name(amp, n.amp)
    ppi = default_name(ppi, n.ppi)

    model.add_reaction(
        name=rxn,
        fn=reversible_michaelis_menten_3s_3p,
        stoichiometry={
            pyruvate: -1,
            atp: -1,
            pep: 1,
        },
        args=[
            pyruvate,
            atp,
            pi,
            pep,
            amp,
            ppi,
            default_vmax(
                model,
                rxn=rxn,
                e0=e0,
                kcat=kcat,
                e0_value=1.0,  # Source
                kcat_value=1.0,  # Source
            ),
            default_kms(model, rxn=rxn, par=kms, value=0.1),
            default_kmp(model, rxn=rxn, par=kmp, value=1.0),
            default_keq(model, rxn=rxn, par=keq, value=0.0096),
        ],
    )
    return model

"""Enzyme building block: tp export."""

from __future__ import annotations

from typing import TYPE_CHECKING

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_name,
    static,
)

if TYPE_CHECKING:
    from mxlpy import Model


def _rate_translocator(
    pi: float,
    pga: float,
    gap: float,
    dhap: float,
    k_pxt: float,
    p_ext: float,
    k_pi: float,
    k_pga: float,
    k_gap: float,
    k_dhap: float,
) -> float:
    return 1 + (1 + k_pxt / p_ext) * (
        pi / k_pi + pga / k_pga + gap / k_gap + dhap / k_dhap
    )


def _rate_out(
    s1: float,
    n_total: float,
    vmax_efflux: float,
    k_efflux: float,
) -> float:
    return vmax_efflux * s1 / (n_total * k_efflux)


def add_pga_exporter(
    model: Model,
    rxn: str,
    pga: str,
    n_translocator: str,
    vmax_export: str,
    km_pga: str,
) -> Model:
    """Add pga exporter to model."""
    rxn = default_name(rxn, n.ex_pga)

    model.add_reaction(
        name=rxn,
        fn=_rate_out,
        stoichiometry={
            pga: -1,
        },
        args=[
            pga,
            n_translocator,
            vmax_export,
            km_pga,
        ],
    )
    return model


def add_gap_exporter(
    model: Model,
    rxn: str,
    gap: str,
    n_translocator: str,
    vmax_export: str,
    km_gap: str,
) -> Model:
    """Add gap exporter to model."""
    rxn = default_name(rxn, n.ex_gap)

    model.add_reaction(
        name=rxn,
        fn=_rate_out,
        stoichiometry={
            gap: -1,
        },
        args=[
            gap,
            n_translocator,
            vmax_export,
            km_gap,
        ],
    )

    return model


def add_dhap_exporter(
    model: Model,
    *,
    rxn: str,
    dhap: str,
    n_translocator: str,
    vmax_export: str,
    km_dhap: str,
) -> Model:
    """Add dhap exporter to model."""
    rxn = default_name(rxn, n.ex_dhap)

    model.add_reaction(
        name=rxn,
        fn=_rate_out,
        stoichiometry={
            dhap: -1,
        },
        args=[
            dhap,
            n_translocator,
            vmax_export,
            km_dhap,
        ],
    )
    return model


def add_triose_phosphate_exporters(
    model: Model,
    *,
    pga_rxn: str | None = None,
    gap_rxn: str | None = None,
    dhap_rxn: str | None = None,
    pi: str | None = None,
    pga: str | None = None,
    gap: str | None = None,
    dhap: str | None = None,
    e0: str | None = None,
    km_pga: str | None = None,
    km_gap: str | None = None,
    km_dhap: str | None = None,
    km_pi_ext: str | None = None,
    km_pi: str | None = None,
    kcat_export: str | None = None,
) -> Model:
    """Add triose phosphate exporters to model."""
    n_translocator = "N_translocator"
    pga_rxn = default_name(pga_rxn, n.ex_pga)
    gap_rxn = default_name(gap_rxn, n.ex_gap)
    dhap_rxn = default_name(dhap_rxn, n.ex_dhap)
    pi = default_name(pi, n.pi)
    pga = default_name(pga, n.pga)
    gap = default_name(gap, n.gap)
    dhap = default_name(dhap, n.dhap)

    pi_ext = static(model, n.pi_ext(), 0.5)

    km_pga = static(model, n.km(pga_rxn), 0.25) if km_pga is None else km_pga
    km_gap = static(model, n.km(gap_rxn), 0.075) if km_gap is None else km_gap
    km_dhap = static(model, n.km(dhap_rxn), 0.077) if km_dhap is None else km_dhap
    km_pi_ext = static(model, n.km(n_translocator, n.pi_ext()), 0.74)
    km_pi = static(model, n.km(n_translocator, n.pi()), 0.63)

    kcat_export = (
        static(model, n.kcat(n_translocator), 0.25 * 8)
        if kcat_export is None
        else kcat_export
    )

    e0 = static(model, n.e0(n_translocator), 1.0) if e0 is None else e0
    model.add_derived(
        vmax_export := n.vmax(pga_rxn), fn=mass_action_1s, args=[kcat_export, e0]
    )

    model.add_derived(
        name=n_translocator,
        fn=_rate_translocator,
        args=[
            n.pi(),
            n.pga(),
            n.gap(),
            n.dhap(),
            km_pi_ext,
            pi_ext,
            km_pi,
            km_pga,
            km_gap,
            km_dhap,
        ],
    )
    add_pga_exporter(
        model=model,
        rxn=pga_rxn,
        pga=pga,
        n_translocator=n_translocator,
        vmax_export=vmax_export,
        km_pga=km_pga,
    )
    add_gap_exporter(
        model=model,
        rxn=gap_rxn,
        gap=gap,
        n_translocator=n_translocator,
        vmax_export=vmax_export,
        km_gap=km_gap,
    )
    add_dhap_exporter(
        model=model,
        rxn=dhap_rxn,
        dhap=dhap,
        n_translocator=n_translocator,
        vmax_export=vmax_export,
        km_dhap=km_dhap,
    )

    return model

"""name.

EC FIXME

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s
from mxlbricks.utils import (
    default_name,
    static,
)


def _rate_state_transition_ps1_ps2(
    ant: float,
    pox: float,
    p_tot: float,
    k_stt7: float,
    km_st: float,
    n_st: float,
) -> float:
    """STT7-kinase phosphorylation of LHC; inhibited by oxidised PQ (state 1 → 2 transition)."""
    return k_stt7 * (1 / (1 + (pox / p_tot / km_st) ** n_st)) * ant


def add_state_transition_12(
    model: Model,
    *,
    rxn: str | None = None,
    lhc: str | None = None,
    pq_ox: str | None = None,
    total_pq: str | None = None,
    kstt7: str | None = None,
    kms: str | None = None,
    n_st: str | None = None,
) -> Model:
    """Add state transition 12 (name) to model."""
    rxn = default_name(rxn, n.lhc_state_transition_12)
    lhc = default_name(lhc, n.lhc)
    pq_ox = default_name(pq_ox, n.pq_ox)
    total_pq = default_name(total_pq, n.total_pq)

    model.add_reaction(
        name=rxn,
        fn=_rate_state_transition_ps1_ps2,
        stoichiometry={
            lhc: -1,
        },
        args=[
            lhc,
            pq_ox,
            total_pq,
            static(model, "kStt7", 0.0035) if kstt7 is None else kstt7,
            static(model, n.km(rxn), 0.2) if kms is None else kms,
            static(model, "n_ST", 2.0) if n_st is None else n_st,
        ],
    )
    return model


def add_state_transition_21(
    model: Model,
    *,
    rxn: str | None = None,
    lhc: str | None = None,
    lhcp: str | None = None,
    kpph: str | None = None,
) -> Model:
    """Add state transition 21 (name) to model."""
    rxn = n.lhc_state_transition_21()
    lhc = default_name(lhc, n.lhc)
    lhcp = default_name(lhcp, n.lhcp)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry={
            lhc: 1,
        },
        args=[
            lhcp,
            static(model, "kPph1", 0.0013) if kpph is None else kpph,
        ],
    )
    return model


def add_state_transitions(
    model: Model,
    *,
    rxn_12: str | None = None,
    rxn_21: str | None = None,
    lhc: str | None = None,
    lhcp: str | None = None,
    pq_ox: str | None = None,
    total_pq: str | None = None,
    kstt7: str | None = None,
    kms: str | None = None,
    n_st: str | None = None,
    kpph: str | None = None,
) -> Model:
    """Add state transitions (name) to model."""
    add_state_transition_12(
        model=model,
        rxn=rxn_12,
        lhc=lhc,
        pq_ox=pq_ox,
        total_pq=total_pq,
        kstt7=kstt7,
        kms=kms,
        n_st=n_st,
    )
    add_state_transition_21(
        model=model,
        rxn=rxn_21,
        lhc=lhc,
        lhcp=lhcp,
        kpph=kpph,
    )

    return model

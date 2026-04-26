"""name.

EC FIXME

Equilibrator
"""

from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import (
    mass_action_1s,
    mass_action_2s,
    michaelis_menten_1s,
    proportional,
)
from mxlbricks.utils import (
    default_name,
    default_par,
    filter_stoichiometry,
)


def add_cbb_pfd_linear_speedup(
    model: Model,
    *,
    der: str | None = None,
    pfd: str | None = None,
    factor: str | None = None,
) -> Model:
    """Add speed-up of CBB enzymes using a linear function."""
    der = default_name(der, n.light_speedup)
    pfd = default_name(pfd, n.pfd)

    model.add_derived(
        der,
        fn=proportional,
        args=[
            pfd,
            default_par(model, par=factor, name=n.kf(der), value=2.0),
        ],
    )
    return model


def add_cbb_pfd_mm_speedup(
    model: Model,
    *,
    der: str | None = None,
    pfd: str | None = None,
    km: str | None = None,
    vmax: str | None = None,
) -> Model:
    """Add speed-up of CBB enzymes using a michaelis-menten curve."""
    der = default_name(der, n.light_speedup)
    pfd = default_name(pfd, n.pfd)

    model.add_derived(
        der,
        fn=michaelis_menten_1s,
        args=[
            pfd,
            default_par(model, par=vmax, name=n.vmax(der), value=6.0),
            default_par(model, par=km, name=n.km(der), value=150.0),
        ],
    )
    return model


def add_fd_tr_reductase_2021(
    model: Model,
    *,
    rxn: str | None = None,
    tr_ox: str | None = None,
    fd_red: str | None = None,
    tr_red: str | None = None,
    fd_ox: str | None = None,
    kf: str | None = None,
) -> Model:
    """Equilibrator.

    Thioredoxin(ox)(aq) + 2 ferredoxin(red)(aq) ⇌ Thioredoxin(red)(aq) + 2 ferredoxin(ox)(aq)
    Keq = 4.9e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25).
    """
    rxn = default_name(rxn, n.ferredoxin_thioredoxin_reductase)
    tr_ox = default_name(tr_ox, n.tr_ox)
    fd_red = default_name(fd_red, n.fd_red)
    tr_red = default_name(tr_red, n.tr_red)
    fd_ox = default_name(fd_ox, n.fd_ox)

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                tr_ox: -1,
                fd_red: -1,
                tr_red: 1,
                fd_ox: 1,
            },
        ),
        args=[
            tr_ox,
            fd_red,
            default_par(model, par=kf, name=n.kf(rxn), value=1.0),
        ],
    )
    return model


def add_fd_tr_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    tr_ox: str | None = None,
    fd_red: str | None = None,
    tr_red: str | None = None,
    fd_ox: str | None = None,
    kf: str | None = None,
) -> Model:
    """Equilibrator.

    Thioredoxin(ox)(aq) + 2 ferredoxin(red)(aq) ⇌ Thioredoxin(red)(aq) + 2 ferredoxin(ox)(aq)
    Keq = 4.9e3 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25).
    """
    rxn = default_name(rxn, n.ferredoxin_thioredoxin_reductase)
    tr_ox = default_name(tr_ox, n.tr_ox)
    fd_red = default_name(fd_red, n.fd_red)
    tr_red = default_name(tr_red, n.tr_red)
    fd_ox = default_name(fd_ox, n.fd_ox)

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                tr_ox: -1,
                fd_red: -2,
                tr_red: 1,
                fd_ox: 2,
            },
        ),
        args=[
            tr_ox,
            fd_red,
            default_par(model, par=kf, name=n.kf(rxn), value=1.0),
        ],
    )
    return model


def add_nadph_tr_reductase(
    model: Model,
    *,
    rxn: str | None = None,
    tr_ox: str | None = None,
    nadph: str | None = None,
    tr_red: str | None = None,
    nadp: str | None = None,
    kf: str | None = None,
) -> Model:
    """Equilibrator.

    Thioredoxin(ox)(aq) + NADPH(aq) ⇌ Thioredoxin(red)(aq) + NADP(aq)
    Keq = 2e1 (@ pH = 7.5, pMg = 3.0, Ionic strength = 0.25).
    """
    rxn = default_name(rxn, n.ferredoxin_thioredoxin_reductase)
    tr_ox = default_name(tr_ox, n.tr_ox)
    nadph = default_name(nadph, n.nadph)
    tr_red = default_name(tr_red, n.tr_red)
    nadp = default_name(nadp, n.nadp)

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                tr_ox: -1,
                nadph: -1,
                tr_red: 1,
                nadp: 1,
            },
        ),
        args=[
            tr_ox,
            nadph,
            default_par(model, par=kf, name=n.kf(rxn), value=1.0),
        ],
    )
    return model


def add_tr_e_activation(
    model: Model,
    *,
    rxn: str | None = None,
    e_inactive: str | None = None,
    tr_red: str | None = None,
    e_active: str | None = None,
    tr_ox: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add tr e activation (name) to model."""
    rxn = default_name(rxn, n.tr_activation)
    e_inactive = default_name(e_inactive, n.e_inactive)
    tr_red = default_name(tr_red, n.tr_red)
    e_active = default_name(e_active, n.e_active)
    tr_ox = default_name(tr_ox, n.tr_ox)

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                e_inactive: -1,
                tr_red: -1,
                e_active: 1,
                tr_ox: 1,
            },
        ),
        args=[
            e_inactive,
            tr_red,
            default_par(model, par=kf, name=n.kf(rxn), value=1.0),
        ],
    )
    return model


def add_tr_e_activation2021(
    model: Model,
    *,
    rxn: str | None = None,
    e_inactive: str | None = None,
    tr_red: str | None = None,
    e_active: str | None = None,
    tr_ox: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add tr e activation2021 (name) to model."""
    rxn = default_name(rxn, n.tr_activation)
    e_inactive = default_name(e_inactive, n.e_inactive)
    tr_red = default_name(tr_red, n.tr_red)
    e_active = default_name(e_active, n.e_active)
    tr_ox = default_name(tr_ox, n.tr_ox)

    model.add_reaction(
        name=rxn,
        fn=mass_action_2s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                e_inactive: -5,
                tr_red: -5,
                e_active: 5,
                tr_ox: 5,
            },
        ),
        args=[
            e_inactive,
            tr_red,
            default_par(model, par=kf, name=n.kf(rxn), value=1.0),
        ],
    )
    return model


def add_e_relaxation(
    model: Model,
    *,
    rxn: str | None = None,
    e_active: str | None = None,
    e_inactive: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add e relaxation (name) to model."""
    rxn = default_name(rxn, n.tr_inactivation)
    e_active = default_name(e_active, n.e_active)
    e_inactive = default_name(e_inactive, n.e_inactive)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                e_active: -1,
                e_inactive: 1,
            },
        ),
        args=[
            e_active,
            default_par(model, par=kf, name=n.kf(rxn), value=0.1),
        ],
    )
    return model


def add_e_relaxation_2021(
    model: Model,
    *,
    rxn: str | None = None,
    e_active: str | None = None,
    e_inactive: str | None = None,
    kf: str | None = None,
) -> Model:
    """Add e relaxation 2021 (name) to model."""
    rxn = default_name(rxn, n.tr_inactivation)
    e_active = default_name(e_active, n.e_active)
    e_inactive = default_name(e_inactive, n.e_inactive)

    model.add_reaction(
        name=rxn,
        fn=mass_action_1s,
        stoichiometry=filter_stoichiometry(
            model,
            {
                e_active: -5,
                e_inactive: 5,
            },
        ),
        args=[
            e_active,
            default_par(model, par=kf, name=n.kf(rxn), value=0.1),
        ],
    )
    return model


def add_thioredoxin_regulation(model: Model) -> Model:
    """Add thioredoxin regulation (name) to model."""
    add_fd_tr_reductase(model)
    add_tr_e_activation(model)
    add_e_relaxation(model)
    return model


def add_thioredoxin_regulation2021(model: Model) -> Model:
    """Add thioredoxin regulation2021 (name) to model."""
    add_fd_tr_reductase_2021(model)
    add_tr_e_activation2021(model)
    add_e_relaxation_2021(model)
    return model


def add_thioredoxin_regulation_from_nadph(model: Model) -> Model:
    """Add thioredoxin regulation from nadph (name) to model."""
    add_nadph_tr_reductase(model)
    add_tr_e_activation(model)
    add_e_relaxation(model)
    return model

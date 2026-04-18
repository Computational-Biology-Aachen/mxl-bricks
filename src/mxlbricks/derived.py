import numpy as np
from mxlpy import Model

from mxlbricks import names as n
from mxlbricks.fns import div, mass_action_1s, moiety_1, moiety_2
from mxlbricks.utils import default_name, default_par


def add_ascorbate_moiety(
    model: Model,
    *,
    name: str | None = None,
    mda: str | None = None,
    dha: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.ascorbate),
        fn=moiety_2,
        args=[
            default_name(mda, n.mda),
            default_name(dha, n.dha),
            default_par(model, par=total, name=n.total_ascorbate(), value=10),
        ],
    )
    return model


def add_adenosin_moiety(
    model: Model,
    *,
    name: str | None = None,
    atp: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.adp),
        fn=moiety_1,
        args=[
            default_name(atp, n.atp),
            default_par(model, par=total, name=n.total_adenosines(), value=0.5),
        ],
    )
    return model


def add_enzyme_moiety(
    model: Model,
    *,
    name: str | None = None,
    e_inactive: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.e_active),
        fn=moiety_1,
        args=[
            default_name(e_inactive, n.e_inactive),
            default_par(model, par=total, name=n.e_total(), value=6.0),
        ],
    )
    return model


def add_ferredoxin_moiety(
    model: Model,
    *,
    name: str | None = None,
    fd_ox: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.fd_red),
        fn=moiety_1,
        args=[
            default_name(fd_ox, n.fd_ox),
            default_par(model, par=total, name=n.total_ferredoxin(), value=5.0),
        ],
    )
    return model


def add_glutamate_moiety(
    model: Model,
    *,
    name: str | None = None,
    glutamate: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.oxoglutarate),
        fn=moiety_1,
        args=[
            default_name(glutamate, n.glutamate),
            default_par(model, par=total, name=n.total_glutamate(), value=3.0),
        ],
    )
    return model


def _glutathion_moiety(
    gssg: float,
    gs_total: float,
) -> float:
    return gs_total - 2 * gssg


def add_glutathion_moiety(
    model: Model,
    *,
    name: str | None = None,
    glutathion_ox: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.glutathion_red),
        fn=_glutathion_moiety,
        args=[
            default_name(glutathion_ox, n.glutathion_ox),
            default_par(model, par=total, name=n.total_glutathion(), value=10.0),
        ],
    )
    return model


def add_hco3_from_co2(
    model: Model,
    *,
    name: str | None = None,
    co2: str | None = None,
    factor: str | None = None,
) -> Model:
    return model.add_derived(
        name=default_name(name, n.hco3),
        fn=mass_action_1s,
        args=[
            default_name(co2, n.co2),
            default_par(model, par=factor, name="CO2/HCO3 ratio", value=50),
        ],
    )


def add_lhc_moiety(
    model: Model,
    *,
    name: str | None = None,
    lhc: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.lhcp),
        fn=moiety_1,
        args=[
            default_name(lhc, n.lhc),
            default_par(model, par=total, name=n.total_lhc(), value=1.0),
        ],
    )
    return model


def add_nad_moiety(
    model: Model,
    *,
    name: str | None = None,
    nadh: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.nad),
        fn=moiety_1,
        args=[
            default_name(nadh, n.nadh),
            default_par(model, par=total, name=n.total_nad(), value=0.86),
        ],
    )
    return model


def add_nadp_moiety(
    model: Model,
    *,
    name: str | None = None,
    nadph: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.nadp),
        fn=moiety_1,
        args=[
            default_name(nadph, n.nadph),
            default_par(model, par=total, name=n.total_nadp(), value=0.5),
        ],
    )
    return model


def add_plastocyanin_moiety(
    model: Model,
    *,
    name: str | None = None,
    pc_ox: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.pc_red),
        fn=moiety_1,
        args=[
            default_name(pc_ox, n.pc_ox),
            default_par(model, par=total, name=n.total_pc(), value=4.0),
        ],
    )
    return model


def add_plastoquinone_moiety(
    model: Model,
    *,
    name: str | None = None,
    pq_ox: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.pq_red),
        fn=moiety_1,
        args=[
            default_name(pq_ox, n.pq_ox),
            default_par(model, par=total, name=n.total_pq(), value=17.5),
        ],
    )
    return model


def add_carotenoid_moiety(
    model: Model,
    *,
    name: str | None = None,
    vx: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.zx),
        fn=moiety_1,
        args=[
            default_name(vx, n.vx),
            default_par(model, par=total, name=n.total_carotenoids(), value=1.0),
        ],
    )
    return model


def add_thioredoxin_moiety(
    model: Model,
    *,
    name: str | None = None,
    tr_ox: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.tr_red),
        fn=moiety_1,
        args=[
            default_name(tr_ox, n.tr_ox),
            default_par(model, par=total, name=n.total_thioredoxin(), value=1.0),
        ],
    )
    return model


def add_psbs_moietry(
    model: Model,
    *,
    name: str | None = None,
    psbs_de: str | None = None,
    total: str | None = None,
) -> Model:
    """Derive protonated form from deprotonated form"""
    model.add_derived(
        name=default_name(name, n.psbs_pr),
        fn=moiety_1,
        args=[
            default_name(psbs_de, n.psbs_de),
            default_par(model, par=total, name=n.total_psbs(), value=1.0),
        ],
    )
    return model


def add_rt(
    model: Model,
    name: str | None = None,
    r: str | None = None,
    t: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(name, n.rt),
        fn=mass_action_1s,
        args=[
            default_par(model, par=r, name="R", value=0.0083),
            default_par(model, par=t, name="T", value=298.0),
        ],
    )
    return model


def _keq_pq_red(
    E0_QA: float,
    F: float,
    E0_PQ: float,
    pHstroma: float,
    dG_pH: float,
    RT: float,
) -> float:
    dg1 = -E0_QA * F
    dg2 = -2 * E0_PQ * F
    dg = -2 * dg1 + dg2 + 2 * pHstroma * dG_pH
    return np.exp(-dg / RT)


def add_plastoquinone_keq(
    model: Model,
    *,
    pq_red: str | None = None,
    ph: str | None = None,
) -> Model:
    model.add_parameter("E^0_QA", -0.14)
    model.add_parameter("E^0_PQ", 0.354)

    model.add_derived(
        n.keq(default_name(pq_red, n.pq_red)),
        _keq_pq_red,
        args=[
            "E^0_QA",
            "F",
            "E^0_PQ",
            default_name(ph, n.ph),
            "dG_pH",
            "RT",
        ],
    )
    return model


def _quencher(
    Psbs: float,
    Vx: float,
    Psbsp: float,
    Zx: float,
    y0: float,
    y1: float,
    y2: float,
    y3: float,
    kZSat: float,
) -> float:
    """co-operative 4-state quenching mechanism
    gamma0: slow quenching of (Vx - protonation)
    gamma1: fast quenching (Vx + protonation)
    gamma2: fastest possible quenching (Zx + protonation)
    gamma3: slow quenching of Zx present (Zx - protonation)
    """
    ZAnt = Zx / (Zx + kZSat)
    return y0 * Vx * Psbs + y1 * Vx * Psbsp + y2 * ZAnt * Psbsp + y3 * ZAnt * Psbs


def add_quencher(
    model: Model,
    *,
    quencher: str | None = None,
    psbs_de: str | None = None,
    vx: str | None = None,
    psbs_pr: str | None = None,
    zx: str | None = None,
) -> Model:
    model.add_parameter("gamma0", 0.1)
    model.add_parameter("gamma1", 0.25)
    model.add_parameter("gamma2", 0.6)
    model.add_parameter("gamma3", 0.15)
    model.add_parameter("kZSat", 0.12)
    model.add_derived(
        name=default_name(quencher, n.quencher),
        fn=_quencher,
        args=[
            default_name(psbs_de, n.psbs_de),
            default_name(vx, n.vx),
            default_name(psbs_pr, n.psbs_pr),
            default_name(zx, n.zx),
            "gamma0",
            "gamma1",
            "gamma2",
            "gamma3",
            "kZSat",
        ],
    )
    return model


def _ph_lumen(protons: float) -> float:
    return -np.log10(protons * 0.00025)


def _dg_ph(r: float, t: float) -> float:
    return np.log(10) * r * t


def add_ph_lumen(
    model: Model,
    *,
    ph: str | None = None,
    r: str | None = None,
    t: str | None = None,
    h: str | None = None,
) -> Model:
    model.add_derived(
        "dG_pH",
        _dg_ph,
        args=[
            default_name(r, lambda: "R"),
            default_name(t, lambda: "T"),
        ],
    )

    model.add_derived(
        name=default_name(ph, lambda: n.ph("_lumen")),
        fn=_ph_lumen,
        args=[
            default_name(h, lambda: n.h("_lumen")),
        ],
    )
    return model


def _pi_cbb(
    phosphate_total: float,
    pga: float,
    bpga: float,
    gap: float,
    dhap: float,
    fbp: float,
    f6p: float,
    g6p: float,
    g1p: float,
    sbp: float,
    s7p: float,
    e4p: float,
    x5p: float,
    r5p: float,
    rubp: float,
    ru5p: float,
    atp: float,
) -> float:
    return phosphate_total - (
        pga
        + 2 * bpga
        + gap
        + dhap
        + 2 * fbp
        + f6p
        + g6p
        + g1p
        + 2 * sbp
        + s7p
        + e4p
        + x5p
        + r5p
        + 2 * rubp
        + ru5p
        + atp
    )


def _pi_cbb_pr(
    phosphate_total: float,
    pga: float,
    bpga: float,
    gap: float,
    dhap: float,
    fbp: float,
    f6p: float,
    g6p: float,
    g1p: float,
    sbp: float,
    s7p: float,
    e4p: float,
    x5p: float,
    r5p: float,
    rubp: float,
    ru5p: float,
    atp: float,
    pgo: float,
) -> float:
    return phosphate_total - (
        pga
        + 2 * bpga
        + gap
        + dhap
        + 2 * fbp
        + f6p
        + g6p
        + g1p
        + 2 * sbp
        + s7p
        + e4p
        + x5p
        + r5p
        + 2 * rubp
        + ru5p
        + atp
        + pgo
    )


def add_orthophosphate_moiety_cbb(
    model: Model,
    *,
    pi: str | None = None,
    total_pi: str | None = None,
    pga: str | None = None,
    bpga: str | None = None,
    gap: str | None = None,
    dhap: str | None = None,
    fbp: str | None = None,
    f6p: str | None = None,
    g6p: str | None = None,
    g1p: str | None = None,
    sbp: str | None = None,
    s7p: str | None = None,
    e4p: str | None = None,
    x5p: str | None = None,
    r5p: str | None = None,
    rubp: str | None = None,
    ru5p: str | None = None,
    atp: str | None = None,
    total: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(pi, n.pi),
        fn=_pi_cbb,
        args=[
            default_par(
                model,
                par=total,
                name=default_name(total_pi, n.total_orthophosphate),
                value=15.0,
            ),
            default_name(pga, n.pga),
            default_name(bpga, n.bpga),
            default_name(gap, n.gap),
            default_name(dhap, n.dhap),
            default_name(fbp, n.fbp),
            default_name(f6p, n.f6p),
            default_name(g6p, n.g6p),
            default_name(g1p, n.g1p),
            default_name(sbp, n.sbp),
            default_name(s7p, n.s7p),
            default_name(e4p, n.e4p),
            default_name(x5p, n.x5p),
            default_name(r5p, n.r5p),
            default_name(rubp, n.rubp),
            default_name(ru5p, n.ru5p),
            default_name(atp, n.atp),
        ],
    )

    return model


def add_orthophosphate_moiety_cbb_pr(
    model: Model,
    *,
    pi: str | None = None,
    total_pi: str | None = None,
    pga: str | None = None,
    bpga: str | None = None,
    gap: str | None = None,
    dhap: str | None = None,
    fbp: str | None = None,
    f6p: str | None = None,
    g6p: str | None = None,
    g1p: str | None = None,
    sbp: str | None = None,
    s7p: str | None = None,
    e4p: str | None = None,
    x5p: str | None = None,
    r5p: str | None = None,
    rubp: str | None = None,
    ru5p: str | None = None,
    atp: str | None = None,
    total: str | None = None,
    pgo: str | None = None,
) -> Model:
    model.add_derived(
        name=default_name(pi, n.pi),
        fn=_pi_cbb,
        args=[
            default_par(
                model,
                par=total,
                name=default_name(total_pi, n.total_orthophosphate),
                value=20.0,
            ),
            default_name(pga, n.pga),
            default_name(bpga, n.bpga),
            default_name(gap, n.gap),
            default_name(dhap, n.dhap),
            default_name(fbp, n.fbp),
            default_name(f6p, n.f6p),
            default_name(g6p, n.g6p),
            default_name(g1p, n.g1p),
            default_name(sbp, n.sbp),
            default_name(s7p, n.s7p),
            default_name(e4p, n.e4p),
            default_name(x5p, n.x5p),
            default_name(r5p, n.r5p),
            default_name(rubp, n.rubp),
            default_name(ru5p, n.ru5p),
            default_name(atp, n.atp),
            default_name(pgo, n.pgo),
        ],
    )

    return model


def _rate_fluorescence(
    Q: float,
    B0: float,
    B2: float,
    ps2cs: float,
    k2: float,
    kF: float,
    kH: float,
) -> float:
    return ps2cs * kF * B0 / (kF + k2 + kH * Q) + ps2cs * kF * B2 / (kF + kH * Q)


def add_readouts(
    model: Model,
    *,
    pq: bool = False,
    fd: bool = False,
    pc: bool = False,
    nadph: bool = False,
    atp: bool = False,
    fluorescence: bool = False,
) -> Model:
    if pq:
        model.add_readout(
            name="PQ_ox/tot",
            fn=div,
            args=[n.pq_red(), n.total_pq()],
        )
    if fd:
        model.add_readout(
            name="Fd_ox/tot",
            fn=div,
            args=[n.fd_red(), n.total_ferredoxin()],
        )
    if pc:
        model.add_readout(
            name="PC_ox/tot",
            fn=div,
            args=[n.pc_red(), n.total_pc()],
        )
    if nadph:
        model.add_readout(
            name="NADPH/tot",
            fn=div,
            args=[n.nadph(), n.total_nadp()],
        )
    if atp:
        model.add_readout(
            name="ATP/tot",
            fn=div,
            args=[n.atp(), n.total_adenosines()],
        )
    if fluorescence:
        model.add_readout(
            name=n.fluorescence(),
            fn=_rate_fluorescence,
            args=[
                n.quencher(),
                n.b0(),
                n.b2(),
                n.ps2cs(),
                "k2",
                "kF",
                "kH",
            ],
        )
    return model

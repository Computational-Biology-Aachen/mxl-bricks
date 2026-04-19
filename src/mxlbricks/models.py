from __future__ import annotations

from typing import Literal

from mxlpy import Model, Parameter, Variable, units
from sympy.physics.units.quantities import Quantity

from . import names as n
from ._tmp import get_matuszynska2016_npq
from .derived import (
    add_adenosin_moiety,
    add_ascorbate_moiety,
    add_carotenoid_moiety,
    add_enzyme_moiety,
    add_ferredoxin_moiety,
    add_glutathion_moiety,
    add_lhc_moiety,
    add_nadp_moiety,
    add_orthophosphate_moiety_cbb,
    add_ph_lumen,
    add_plastocyanin_moiety,
    add_plastoquinone_keq,
    add_plastoquinone_moiety,
    add_psbs_moietry,
    add_quencher,
    add_readouts,
    add_rt,
    add_thioredoxin_moiety,
)
from .enzymes import (
    add_aldolase_dhap_e4p_req,
    add_aldolase_dhap_gap_req,
    add_ascorbate_peroxidase,
    add_atp_consumption,
    add_atp_synthase_mm,
    add_atp_synthase_mmol_chl,
    add_atp_synthase_static_protons,
    add_b6f,
    add_catalase,
    add_cbb_pfd_linear_speedup,
    add_cbb_pfd_mm_speedup,
    add_cyclic_electron_flow,
    add_dehydroascorbate_reductase,
    add_fbpase,
    add_ferredoxin_reductase,
    add_fnr_mm,
    add_fnr_mmol_chl,
    add_g1p_efflux,
    add_gadph,
    add_glucose_6_phosphate_isomerase_re,
    add_glutathion_reductase_irrev,
    add_glycine_decarboxylase_yokota,
    add_glycine_transaminase_yokota,
    add_glycolate_oxidase_yokota,
    add_hpa_outflux,
    add_lhc_deprotonation,
    add_lhc_protonation,
    add_mda_reductase1,
    add_mda_reductase2,
    add_nadph_consumption,
    add_ndh,
    add_phosphoglucomutase,
    add_phosphoglycerate_kinase_poolman,
    add_phosphoglycolate_influx,
    add_phosphoribulokinase,
    add_photosystems,
    add_proton_leak,
    add_ps2_cross_section,
    add_ptox,
    add_ribose_5_phosphate_isomerase,
    add_ribulose_5_phosphate_3_epimerase,
    add_rubisco_poolman,
    add_sbpase,
    add_serine_glyoxylate_transaminase_irreversible,
    add_state_transitions,
    add_thioredoxin_regulation2021,
    # add_thioredoxin_regulation_from_nadph,
    add_transketolase_x5p_e4p_f6p_gap,
    add_transketolase_x5p_r5p_s7p_gap,
    add_triose_phosphate_exporters,
    add_triose_phosphate_isomerase,
    add_violaxanthin_epoxidase,
    add_zeaxanthin_epoxidase,
)
from .utils import fcbb_regulated, static, thioredixon_regulated

__all__ = [
    "get_matuszynska2016_npq",
    "get_yokota1985",
    "get_poolman2000",
    "get_matuszynska2016_phd",
    "get_matuszynska2019",
    "get_saadat2021",
    "get_cbb_regulated",
]

mol_chl = Quantity("mol_chl", abbrev="mol_chl")
mmol_mol_chl = units.mmol / mol_chl


def get_yokota1985() -> Model:
    model = Model()
    model.add_variables(
        {
            n.glycolate(): 0.09,
            n.glyoxylate(): 0.7964601770483386,
            n.glycine(): 8.999999999424611,
            n.serine(): 2.5385608670239126,
            n.hydroxypyruvate(): 0.009782608695111009,
            n.h2o2(): 0.010880542843616855,
        }
    )

    add_phosphoglycolate_influx(model=model)
    add_glycolate_oxidase_yokota(model=model)
    add_glycine_transaminase_yokota(model=model)
    add_glycine_decarboxylase_yokota(
        model=model,
        e0=static(model, n.e0(n.glycine_decarboxylase()), 0.5),
    )
    add_serine_glyoxylate_transaminase_irreversible(model=model)
    add_hpa_outflux(model=model)
    add_catalase(model=model)
    return model


def get_poolman2000() -> Model:
    model = Model()
    model.add_variables(
        {
            n.pga(): Variable(0.6387788347932627, unit=units.mmol),
            n.bpga(): Variable(0.0013570885908749779, unit=units.mmol),
            n.gap(): Variable(0.011259431827358068, unit=units.mmol),
            n.dhap(): Variable(0.24770748227012374, unit=units.mmol),
            n.fbp(): Variable(0.01980222074817044, unit=units.mmol),
            n.f6p(): Variable(1.093666906864421, unit=units.mmol),
            n.g6p(): Variable(2.5154338857582377, unit=units.mmol),
            n.g1p(): Variable(0.14589516537322303, unit=units.mmol),
            n.sbp(): Variable(0.09132688566151095, unit=units.mmol),
            n.s7p(): Variable(0.23281380022778891, unit=units.mmol),
            n.e4p(): Variable(0.02836065066520614, unit=units.mmol),
            n.x5p(): Variable(0.03647242425941113, unit=units.mmol),
            n.r5p(): Variable(0.06109130988031577, unit=units.mmol),
            n.rubp(): Variable(0.2672164362349537, unit=units.mmol),
            n.ru5p(): Variable(0.0244365238237522, unit=units.mmol),
            n.atp(): Variable(0.43633201706180874, unit=units.mmol),
        }
    )
    model.add_parameters(
        {
            n.co2(): Parameter(0.2, unit=units.mmol),
            n.nadph(): Parameter(0.21, unit=units.mmol),
            n.h(): Parameter(1.2589254117941661e-05, unit=units.mmol),
        }
    )

    add_adenosin_moiety(
        model,
        total=static(model, n.total_adenosines(), 0.5),
    )
    add_nadp_moiety(
        model,
        total=static(model, n.total_nadp(), 0.5),
    )
    add_orthophosphate_moiety_cbb(
        model,
        total=static(model, n.total_orthophosphate(), 15.0),
    )

    # Reactions
    add_rubisco_poolman(model)
    add_phosphoglycerate_kinase_poolman(model)
    add_gadph(model)
    add_triose_phosphate_isomerase(model)
    add_aldolase_dhap_gap_req(model)
    add_aldolase_dhap_e4p_req(model)
    add_fbpase(model)
    add_transketolase_x5p_e4p_f6p_gap(model)
    add_transketolase_x5p_r5p_s7p_gap(model)
    add_sbpase(model)
    add_ribose_5_phosphate_isomerase(model)
    add_ribulose_5_phosphate_3_epimerase(model)
    add_phosphoribulokinase(model)
    add_glucose_6_phosphate_isomerase_re(model)
    add_phosphoglucomutase(model)
    add_triose_phosphate_exporters(model)
    add_g1p_efflux(model)

    # Other
    add_atp_synthase_static_protons(model)
    return model


def get_matuszynska2016_phd(
    *,
    mode: Literal["matrix", "analytical"] = "matrix",
    chl_lumen: str = "_lumen",
) -> Model:
    model = Model()
    model.add_variables(
        {
            n.atp(): 1.6999999999999997,
            n.pq_ox(): 4.706348349506148,
            n.pc_ox(): 3.9414515288091567,
            n.fd_ox(): 3.7761613271207324,
            n.h(chl_lumen): 7.737821100836988,
            n.lhc(): 0.5105293511676007,
            n.psbs_de(): 0.5000000001374878,
            n.vx(): 0.09090909090907397,
        }
    )
    model.add_parameters(
        {
            n.ph(): 7.9,  # stroma ph
            n.pfd(): 100.0,
            n.nadph(): 0.6,
            n.o2(chl_lumen): 8.0,
            "bH": 100.0,
            "F": 96.485,
            "E^0_PC": 0.38,
            "E^0_P700": 0.48,
            "E^0_FA": -0.55,
            "E^0_Fd": -0.43,
            "E^0_NADP": -0.113,
        }
    )

    add_nadp_moiety(model, total=static(model, n.total_nadp(), 0.8))

    # Moieties / derived compounds
    add_rt(model)
    add_adenosin_moiety(
        model,
        total=static(model, n.total_adenosines(), value=2.55, unit=mmol_mol_chl),
    )
    add_ph_lumen(model)
    add_carotenoid_moiety(model)
    add_ferredoxin_moiety(model)
    add_plastocyanin_moiety(model)
    add_psbs_moietry(model)
    add_lhc_moiety(model)
    add_quencher(model)
    add_plastoquinone_keq(model)
    add_plastoquinone_moiety(model)
    add_ps2_cross_section(model)

    # Reactions
    add_atp_synthase_mmol_chl(model, bh="bH")
    add_b6f(model, bh="bH")
    add_lhc_protonation(model)
    add_lhc_deprotonation(model)
    add_cyclic_electron_flow(model)
    add_violaxanthin_epoxidase(model)
    add_zeaxanthin_epoxidase(model)
    add_fnr_mmol_chl(model)
    add_ndh(model)
    add_photosystems(model, mode=mode, mehler=False)
    add_proton_leak(model)
    add_ptox(model)
    add_state_transitions(model)

    # Misc
    add_atp_consumption(
        model,
        kf=static(model, n.kf(n.ex_atp()), 10.0),
    )
    add_readouts(
        model,
        pq=True,
        fd=True,
        pc=True,
        atp=True,
        fluorescence=True,
    )
    return model


def get_matuszynska2019(
    *,
    variant: Literal["linear-speedup", "mm-speedup"] | None = None,
    mode: Literal["matrix", "analytical"] = "matrix",
    chl_lumen: str = "_lumen",
) -> Model:
    model = Model()
    model.add_variables(
        {
            n.pga(): 0.9928653922138561,
            n.bpga(): 0.0005297732935310749,
            n.gap(): 0.0062663539939955834,
            n.dhap(): 0.13785977143668732,
            n.fbp(): 0.006133532145409954,
            n.f6p(): 0.31271973359685457,
            n.g6p(): 0.719255387166192,
            n.g1p(): 0.041716812452951633,
            n.sbp(): 0.013123745088361893,
            n.s7p(): 0.15890073845176905,
            n.e4p(): 0.007322797350442026,
            n.x5p(): 0.022478763225333428,
            n.r5p(): 0.037651927659696716,
            n.rubp(): 0.13184790283048484,
            n.ru5p(): 0.015060770937455408,
            n.atp(): 1.612922506604933,
            n.fd_ox(): 3.8624032084329674,
            n.h(chl_lumen): 0.002208423037307405,
            n.lhc(): 0.80137477470646,
            n.nadph(): 0.491395685599137,
            n.pc_ox(): 1.885391998090184,
            n.pq_ox(): 10.991562708096392,
            n.psbs_de(): 0.9610220887579118,
            n.vx(): 0.9514408605906095,
        }
    )
    model.add_parameters(
        {
            n.h(): 1.2589254117941661e-05,
            n.ph(): 7.9,
            n.co2(): 0.2,
            n.o2(chl_lumen): 8.0,
            n.pfd(): 100.0,
            "bH": 100.0,
            "F": 96.485,
            "E^0_PC": 0.38,
            "E^0_P700": 0.48,
            "E^0_FA": -0.55,
            "E^0_Fd": -0.43,
            "E^0_NADP": -0.113,
            "convf": 3.2e-2,
        }
    )
    if variant == "linear-speedup":
        add_cbb_pfd_linear_speedup(model)
    elif variant == "mm-speedup":
        add_cbb_pfd_mm_speedup(model)

    # Moieties / derived compounds
    add_rt(model)
    add_ph_lumen(model)
    add_carotenoid_moiety(model)
    add_ferredoxin_moiety(model)
    add_plastocyanin_moiety(model)
    add_psbs_moietry(model)
    add_lhc_moiety(model)
    add_quencher(model)
    add_plastoquinone_keq(model)
    add_plastoquinone_moiety(model)
    add_ps2_cross_section(model)
    add_nadp_moiety(
        model,
        total=static(model, n.total_nadp(), 0.8),
    )
    add_adenosin_moiety(
        model,
        total=static(model, n.total_adenosines(), value=2.55, unit=units.mmol),
    )
    add_orthophosphate_moiety_cbb(
        model,
        total=static(model, n.total_orthophosphate(), 17.05),
    )

    # Reactions
    add_atp_synthase_mm(model, bh="bH", convf="convf")
    add_b6f(model, bh="bH")
    add_lhc_protonation(model)
    add_lhc_deprotonation(model)
    add_cyclic_electron_flow(model)
    add_violaxanthin_epoxidase(model)
    add_zeaxanthin_epoxidase(model)
    add_fnr_mm(model, convf="convf")
    add_ndh(model)
    add_photosystems(
        model,
        mode=mode,
        mehler=False,
        convf="convf",
    )
    add_proton_leak(model)
    add_ptox(model)
    add_state_transitions(model)
    add_rubisco_poolman(
        model,
        e0=(
            fcbb_regulated(model, n.e0(n.rubisco()), 1.0)
            if variant is not None
            else static(model, n.e0(n.rubisco()), 1.0)
        ),
    )
    add_phosphoglycerate_kinase_poolman(model)
    add_gadph(model)
    add_triose_phosphate_isomerase(model)
    add_aldolase_dhap_gap_req(model)
    add_aldolase_dhap_e4p_req(model)
    add_fbpase(
        model,
        e0=(
            fcbb_regulated(model, n.e0(n.fbpase()), 1.0)
            if variant is not None
            else static(model, n.e0(n.fbpase()), 1.0)
        ),
    )
    add_transketolase_x5p_e4p_f6p_gap(model)
    add_transketolase_x5p_r5p_s7p_gap(model)
    add_sbpase(
        model,
        e0=(
            fcbb_regulated(model, n.e0(n.sbpase()), 1.0)
            if variant is not None
            else static(model, n.e0(n.sbpase()), 1.0)
        ),
    )
    add_ribose_5_phosphate_isomerase(model)
    add_ribulose_5_phosphate_3_epimerase(model)
    add_phosphoribulokinase(
        model,
        e0=(
            fcbb_regulated(model, n.e0(n.phosphoribulokinase()), 1.0)
            if variant is not None
            else static(model, n.e0(n.phosphoribulokinase()), 1.0)
        ),
    )
    add_glucose_6_phosphate_isomerase_re(model)
    add_phosphoglucomutase(model)
    add_triose_phosphate_exporters(model)
    add_g1p_efflux(
        model,
        e0=(
            fcbb_regulated(model, n.e0(n.ex_g1p()), 1.0)
            if variant is not None
            else static(model, n.e0(n.ex_g1p()), 1.0)
        ),
    )

    add_readouts(
        model,
        pq=True,
        fd=True,
        pc=True,
        nadph=True,
        atp=True,
        fluorescence=True,
    )
    return model


def get_saadat2021(
    *,
    mode: Literal["matrix", "analytical"] = "matrix",
    chl_lumen: str = "_lumen",
) -> Model:
    model = Model()
    model.add_variables(
        {
            n.pga(): 0.9167729479368978,
            n.bpga(): 0.0003814495319659031,
            n.gap(): 0.00580821050261484,
            n.dhap(): 0.1277806166216142,
            n.fbp(): 0.005269452472931973,
            n.f6p(): 0.2874944558066638,
            n.g6p(): 0.6612372482712676,
            n.g1p(): 0.03835176039761378,
            n.sbp(): 0.011101373736607443,
            n.s7p(): 0.1494578301900007,
            n.e4p(): 0.00668295494870102,
            n.x5p(): 0.020988553174809618,
            n.r5p(): 0.035155825913785584,
            n.rubp(): 0.11293260727162346,
            n.ru5p(): 0.014062330254191594,
            n.atp(): 1.4612747767895344,
            n.fd_ox(): 3.715702384326767,
            n.h(chl_lumen): 0.002086128887296243,
            n.lhc(): 0.7805901436176024,
            n.nadph(): 0.5578718406315588,
            n.pc_ox(): 1.8083642974980014,
            n.pq_ox(): 10.251099271612473,
            n.psbs_de(): 0.9667381262477079,
            n.vx(): 0.9629870646993118,
            n.mda(): 2.0353396709300447e-07,
            n.h2o2(): 1.2034405327140102e-07,
            n.dha(): 1.0296456279861962e-11,
            n.glutathion_ox(): 4.99986167652437e-12,
            # New
            n.tr_ox(): 0.9334426859846461,
            n.e_inactive(): 3.6023635680406634,
        }
    )
    model.add_parameters(
        {
            n.pfd(): 100.0,
            n.co2(): 0.2,
            n.o2(chl_lumen): 8.0,
            n.ph(): 7.9,
            n.h(): 1.2589254117941661e-05,
            "bH": 100.0,
            "F": 96.485,
            "E^0_PC": 0.38,
            "E^0_P700": 0.48,
            "E^0_FA": -0.55,
            "E^0_Fd": -0.43,
            "E^0_NADP": -0.113,
            "convf": 3.2e-2,
        }
    )

    # Moieties / derived compounds
    add_rt(model)
    add_ph_lumen(model)
    add_carotenoid_moiety(model)
    add_ferredoxin_moiety(model)
    add_plastocyanin_moiety(model)
    add_psbs_moietry(model)
    add_lhc_moiety(model)
    add_quencher(model)
    add_plastoquinone_keq(model)
    add_plastoquinone_moiety(model)
    add_ps2_cross_section(model)
    add_thioredoxin_moiety(model)
    add_enzyme_moiety(model)
    add_nadp_moiety(
        model,
        total=static(model, n.total_nadp(), 0.8),
    )
    add_adenosin_moiety(
        model,
        total=static(model, n.total_adenosines(), value=2.55, unit=units.mmol),
    )
    add_orthophosphate_moiety_cbb(
        model,
        total=static(model, n.total_orthophosphate(), 17.05),
    )
    add_thioredoxin_regulation2021(model)
    add_ascorbate_moiety(model)
    add_glutathion_moiety(model)

    # Reactions
    ## PETC
    add_atp_synthase_mm(model, convf="convf", bh="bH")
    add_b6f(model, bh="bH")
    add_lhc_protonation(model)
    add_lhc_deprotonation(model)
    add_cyclic_electron_flow(model)
    add_violaxanthin_epoxidase(model)
    add_zeaxanthin_epoxidase(model)
    add_fnr_mm(
        model,
        convf="convf",
    )
    add_ndh(model)
    add_photosystems(
        model,
        mode=mode,
        mehler=True,
        convf="convf",
    )
    add_ferredoxin_reductase(
        model,
        keq=n.keq(n.ferredoxin_reductase()),  # from add_photosystems
    )
    add_proton_leak(model)
    add_ptox(model)
    add_state_transitions(model)

    ## CBB
    add_rubisco_poolman(
        model,
        e0=thioredixon_regulated(model, n.e0(n.rubisco()), 1.0),
    )
    add_phosphoglycerate_kinase_poolman(model)
    add_gadph(model)
    add_triose_phosphate_isomerase(model)
    add_aldolase_dhap_gap_req(model)
    add_aldolase_dhap_e4p_req(model)
    add_fbpase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.fbpase()), 1.0),
    )
    add_transketolase_x5p_e4p_f6p_gap(model)
    add_transketolase_x5p_r5p_s7p_gap(model)
    add_sbpase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.sbpase()), 1.0),
    )
    add_ribose_5_phosphate_isomerase(model)
    add_ribulose_5_phosphate_3_epimerase(model)
    add_phosphoribulokinase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.phosphoribulokinase()), 1.0),
    )
    add_glucose_6_phosphate_isomerase_re(model)
    add_phosphoglucomutase(model)
    add_triose_phosphate_exporters(model)
    add_g1p_efflux(
        model,
        e0=thioredixon_regulated(model, n.e0(n.ex_g1p()), 1.0),
    )

    ## Mehler
    add_mda_reductase1(model)
    add_mda_reductase2(model)
    add_ascorbate_peroxidase(model)
    add_glutathion_reductase_irrev(model)
    add_dehydroascorbate_reductase(model)

    # Misc
    add_atp_consumption(
        model,
        kf=static(model, n.kf(n.ex_atp()), 0.2),
    )
    add_nadph_consumption(
        model,
        kf=static(model, n.kf(n.ex_nadph()), 0.2),
    )
    add_readouts(
        model,
        pq=True,
        fd=True,
        pc=True,
        nadph=True,
        atp=True,
        fluorescence=True,
    )
    return model


def get_cbb_regulated() -> Model:
    model = Model()
    model.add_variables(
        {
            n.pga(): 0.9167729479368978,
            n.bpga(): 0.0003814495319659031,
            n.gap(): 0.00580821050261484,
            n.dhap(): 0.1277806166216142,
            n.fbp(): 0.005269452472931973,
            n.f6p(): 0.2874944558066638,
            n.g6p(): 0.6612372482712676,
            n.g1p(): 0.03835176039761378,
            n.sbp(): 0.011101373736607443,
            n.s7p(): 0.1494578301900007,
            n.e4p(): 0.00668295494870102,
            n.x5p(): 0.020988553174809618,
            n.r5p(): 0.035155825913785584,
            n.rubp(): 0.11293260727162346,
            n.ru5p(): 0.014062330254191594,
            # n.e_active(): 2.3976364319593366,
            # Energy equivalents
            n.atp(): 1.4612747767895344,
            n.nadph(): 0.5578718406315588,
        }
    )
    model.add_parameters(
        {
            n.pfd(): 100.0,
            n.co2(): 0.2,
            n.h(): 1.2589254117941661e-05,
        }
    )

    # Moieties / derived compounds
    add_nadp_moiety(
        model,
        total=static(model, n.total_nadp(), 0.8),
    )
    add_adenosin_moiety(
        model,
        total=static(model, n.total_adenosines(), value=2.55, unit=units.mmol),
    )
    add_orthophosphate_moiety_cbb(
        model,
        total=static(model, n.total_orthophosphate(), 17.05),
    )

    # CBB
    add_rubisco_poolman(
        model,
        e0=thioredixon_regulated(model, n.e0(n.rubisco()), 1.0),
    )
    add_phosphoglycerate_kinase_poolman(model)
    add_gadph(model)
    add_triose_phosphate_isomerase(model)
    add_aldolase_dhap_gap_req(model)
    add_aldolase_dhap_e4p_req(model)
    add_fbpase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.fbpase()), 1.0),
    )
    add_transketolase_x5p_e4p_f6p_gap(model)
    add_transketolase_x5p_r5p_s7p_gap(model)
    add_sbpase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.sbpase()), 1.0),
    )
    add_ribose_5_phosphate_isomerase(model)
    add_ribulose_5_phosphate_3_epimerase(model)
    add_phosphoribulokinase(
        model,
        e0=thioredixon_regulated(model, n.e0(n.phosphoribulokinase()), 1.0),
    )
    add_glucose_6_phosphate_isomerase_re(model)
    add_phosphoglucomutase(model)
    add_triose_phosphate_exporters(model)
    add_g1p_efflux(
        model,
        e0=thioredixon_regulated(model, n.e0(n.ex_g1p()), 1.0),
    )

    # Misc
    add_atp_consumption(
        model,
        kf=static(model, n.kf(n.ex_atp()), 0.2),
    )
    add_nadph_consumption(
        model,
        kf=static(model, n.kf(n.ex_nadph()), 0.2),
    )
    add_readouts(
        model,
        nadph=True,
        atp=True,
    )
    return model

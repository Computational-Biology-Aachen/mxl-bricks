# ## Saadat 2021
# https://doi.org/10.3389/fpls.2021.750580


from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import Assimulo, Simulator, scan

from mxlbricks import names as n
from mxlbricks.models import get_saadat2021

lumen = "_lumen"

CWD = Path(__file__).parent

names = {
    "PQ": n.pq_ox(),
    "PQred": n.pq_red(),
    "PCred": n.pc_red(),
    "PC": n.pc_ox(),
    "Fd": n.fd_ox(),
    "Fdred": n.fd_red(),
    "TR_ox": n.tr_ox(),
    "TR_red": n.tr_red(),
    "ATP": n.atp(),
    "ADP": n.adp(),
    "NADPH": n.nadph(),
    "NADP": n.nadp(),
    "H": n.h(lumen),
    "pH": n.ph(lumen),
    "LHC": n.lhc(),
    "LHCp": n.lhcp(),
    "Psbs": n.psbs_de(),
    "Psbsp": n.psbs_pr(),
    "Vx": n.vx(),
    "Zx": n.zx(),
    "PGA": n.pga(),
    "BPGA": n.bpga(),
    "GAP": n.gap(),
    "DHAP": n.dhap(),
    "FBP": n.fbp(),
    "F6P": n.f6p(),
    "G6P": n.g6p(),
    "G1P": n.g1p(),
    "SBP": n.sbp(),
    "S7P": n.s7p(),
    "E4P": n.e4p(),
    "X5P": n.x5p(),
    "R5P": n.r5p(),
    "RUBP": n.rubp(),
    "RU5P": n.ru5p(),
    "MDA": n.mda(),
    "H2O2": n.h2o2(),
    "DHA": n.dha(),
    "GSSG": n.glutathion_ox(),
    "GSH": n.glutathion_red(),
    "E_active": n.e_active(),
    "E_inactive": n.e_inactive(),
    "ps2cs": n.ps2cs(),
    "Q": n.quencher(),
    "B0": n.b0(),
    "B1": n.b1(),
    "B2": n.b2(),
    "B3": n.b3(),
    "A0": n.a0(),
    "A1": n.a1(),
    "A2": n.a2(),
    "Fluo": n.fluorescence(),
    "Pi": n.pi(),
    "N": "N_translocator",
    "Keq_ATPsynthase": n.keq(n.atp_synthase()),
    "Keq_B6f": n.keq(n.b6f()),
    "ASC": n.ascorbate(),
    "V1": n.vmax(n.rubisco_carboxylase()),
    "V6": n.vmax(n.fbpase()),
    "V9": n.vmax(n.sbpase()),
    "V13": n.vmax(n.phosphoribulokinase()),
    "Vst": n.vmax(n.ex_g1p()),
    "PQ_redoxstate": "PQ_ox/tot",
    "Fd_redoxstate": "Fd_ox/tot",
    "PC_redoxstate": "PC_ox/tot",
    "NADP_redoxstate": "NADPH/tot",
    "ATP_norm": "ATP/tot",
    "vPS2": n.ps2(),
    "vPS1": n.ps1(),
    "vPTOX": n.ptox(),
    "vNDH": n.ndh(),
    "vB6f": n.b6f(),
    "vCyc": n.cyclic_electron_flow(),
    "vFNR": n.fnr(),
    "vLeak": n.proton_leak(),
    "vSt12": n.lhc_state_transition_12(),
    "vSt21": n.lhc_state_transition_21(),
    "vATPsynthase": n.atp_synthase(),
    "vDeepox": n.violaxanthin_deepoxidase(),
    "vEpox": n.zeaxanthin_epoxidase(),
    "vLhcprotonation": n.lhc_protonation(),
    "vLhcdeprotonation": n.lhc_deprotonation(),
    "vRuBisCO": n.rubisco_carboxylase(),
    "vPGA_kinase": n.phosphoglycerate_kinase(),
    "vBPGA_dehydrogenase": n.gadph(),
    "vTPI": n.triose_phosphate_isomerase(),
    "vAldolase": n.aldolase_dhap_gap(),
    "vFBPase": n.fbpase(),
    "vF6P_Transketolase": n.transketolase_gap_f6p(),
    "v8": n.aldolase_dhap_e4p(),
    "v9": n.sbpase(),
    "v10": n.transketolase_gap_s7p(),
    "v11": n.ribose_phosphate_isomerase(),
    "v12": n.ribulose_phosphate_epimerase(),
    "v13": n.phosphoribulokinase(),
    "vG6P_isomerase": n.g6pi(),
    "vPhosphoglucomutase": n.phosphoglucomutase(),
    "vpga": n.ex_pga(),
    "vgap": n.ex_gap(),
    "vdhap": n.ex_dhap(),
    "vStarch": n.ex_g1p(),
    "vFdred": n.ferredoxin_reductase(),
    "vAscorbate": n.ascorbate_peroxidase(),
    "vMDAreduct": n.mda_reductase2(),
    "vMehler": n.mehler(),
    "vGR": n.glutathion_reductase(),
    "vDHAR": n.dehydroascorbate_reductase(),
    "v3ASC": n.mda_reductase1(),
    "vEX_ATP": n.ex_atp(),
    "vEX_NADPH": n.ex_nadph(),
    "vFdTrReductase": n.ferredoxin_thioredoxin_reductase(),
    "vE_activation": n.tr_activation(),
    "vE_inactivation": n.tr_inactivation(),
}


def test_steady_state() -> None:
    model = get_saadat2021()
    res = (
        Simulator(model, integrator=Assimulo).simulate(100).get_result().unwrap_or_err()
    )
    assert res is not None

    pd.testing.assert_series_equal(
        pd.Series(model.get_initial_conditions()),
        pd.Series(res.get_new_y0()),
    )


def test_steady_state_by_pfd() -> None:
    reference = pd.read_csv(CWD / "ss-by-pfd.csv", index_col=0).rename(columns=names)
    reference.index.name = "PPFD"

    res = scan.steady_state(
        get_saadat2021(), to_scan=pd.DataFrame({n.pfd(): np.arange(100, 1500, 100)})
    ).get_args(
        include_variables=True,
        include_parameters=True,
        include_derived_parameters=True,
        include_derived_variables=True,
        include_readouts=True,
        include_surrogate_fluxes=True,
        include_surrogate_variables=True,
    )

    # Check if anything isn't found
    assert len(reference.columns.difference(res.columns)) == 0

    # Check if all mapped are the same
    to_compare = list(reference.columns.intersection(res.columns))
    pd.testing.assert_frame_equal(reference, res.loc[:, to_compare], rtol=1e-4)

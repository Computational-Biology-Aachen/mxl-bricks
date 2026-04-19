# Matuszyńska 2019
# https://doi.org/10.1111/ppl.12962


from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from mxlpy import Scipy, scan

from mxlbricks import get_matuszynska2019 as get_model
from mxlbricks import names as n

CWD = Path(__file__).parent
lumen = "_lumen"

names = {
    "PQ": n.pq_ox(),
    "PQred": n.pq_red(),
    "PCred": n.pc_red(),
    "PC": n.pc_ox(),
    "Fd": n.fd_ox(),
    "Fdred": n.fd_red(),
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
    "ps2cs": n.ps2cs(),
    "Q": n.quencher(),
    "B0": n.b0(),
    "B1": n.b1(),
    "B2": n.b2(),
    "B3": n.b3(),
    "A1": n.a1(),
    "Fluo": n.fluorescence(),
    "Pi": n.pi(),
    "N": "N_translocator",
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
}


def test_steady_state_by_pfd() -> None:
    reference = pd.read_csv(CWD / "ss-by-pfd.csv", index_col=0).rename(columns=names)
    reference.index.name = "PPFD"

    res = scan.steady_state(
        get_model(variant=None),
        to_scan=pd.DataFrame({n.pfd(): np.arange(100, 1400, 100)}),
        integrator=partial(Scipy, method="Radau"),
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
    pd.testing.assert_frame_equal(reference, res.loc[:, to_compare], rtol=1e-2)


@pytest.mark.extensive
def test_steady_state_by_pfd_analytical() -> None:
    reference = pd.read_csv(CWD / "ss-by-pfd.csv", index_col=0).rename(columns=names)
    reference.index.name = "PPFD"

    res = scan.steady_state(
        get_model(variant=None, mode="matrix"),
        to_scan=pd.DataFrame({n.pfd(): np.arange(100, 1400, 100)}),
        integrator=partial(Scipy, method="Radau"),
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
    pd.testing.assert_frame_equal(reference, res.loc[:, to_compare], rtol=1e-2)

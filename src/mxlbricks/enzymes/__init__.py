"""Enzyme reaction building blocks for mxlpy models."""

from .a5p_isomerase import add_a5p_isomerase
from .acetoacetate_coa_ligase import add_acetoacetate_coa_ligase
from .acetyl_coa_acetyltransfer import add_acetyl_coa_acetyltransfer
from .acetyl_coa_carboxyltransfer import (
    add_acetyl_coa_carboxyltransfer,
    add_acetyl_coa_carboxyltransfer_1i,
)
from .aldolase_dhap_a5p_obp import add_aldolase_dhap_a5p_req
from .aldolase_dhap_e4p_sbp import add_aldolase_dhap_e4p_req
from .aldolase_dhap_gap_fbp import add_aldolase_dhap_gap_req
from .aldolase_gah_dhap_r1p import add_r1p_aldolase
from .aldolase_gah_gap_a5p import add_a5p_aldolase
from .aldolase_glyoyx_gly_ha import add_hydroxyaspartate_aldolase
from .ascorbate_perixodase import add_ascorbate_peroxidase
from .aspartate_aminotransferase import add_aspartate_aminotransferase
from .aspartate_nad_oxidoreductase import add_aspartate_oxidoreductase
from .aspartate_nadp_oxidoreductase import add_aspartate_nadp_oxidoreductase
from .atp_synthase import (
    add_atp_synthase_2024,
    add_atp_synthase_energy_dependent,
    add_atp_synthase_mm,
    add_atp_synthase_mmol_chl,
    add_atp_synthase_static_protons,
)
from .b6f import add_b6f, add_b6f_2024
from .carbonic_anhydrase import add_carbonic_anhydrase_mass_action
from .catalase import add_catalase
from .cef import add_cyclic_electron_flow
from .co2_dissolving import add_co2_dissolving
from .dehydroascorbate_reductase import add_dehydroascorbate_reductase
from .e4p_epimerase import add_e4p_epimerase
from .e4p_isomerase import add_e4p_isomerase
from .enolase import add_enolase
from .erythrolose_kinase import add_erythrulose_kinase
from .ex_atp import add_atp_consumption
from .ex_g1p import add_g1p_efflux
from .ex_nadph import add_nadph_consumption
from .fbpase import add_fbpase
from .fd_reductase import add_ferredoxin_reductase
from .fnr import add_fnr_energy_dependent, add_fnr_mm, add_fnr_mmol_chl, add_fnr_static
from .formate_dehydrogenase import add_formate_dehydrogenase
from .formate_thf_ligase import add_formate_thf_ligase
from .g6pi import add_glucose_6_phosphate_isomerase_re
from .gadph import add_gadph
from .glutamate_dehydrogenase import add_glutamate_dehydrogenase
from .glutathion_reductase import add_glutathion_reductase_irrev
from .glycerate_dehydrogenase import add_glycerate_dehydrogenase, add_hpa_outflux
from .glycerate_kinase import add_glycerate_kinase
from .glycine_decarboxylase import (
    add_glycine_decarboxylase,
    add_glycine_decarboxylase_irreversible,
    add_glycine_decarboxylase_yokota,
)
from .glycine_transaminase import (
    add_glycine_transaminase,
    add_glycine_transaminase_irreversible,
    add_glycine_transaminase_yokota,
)
from .glycoaldehyde_dehydrogenase import add_glycolaldehyde_dehydrogenase
from .glycolate_dehydrogenase import add_glycolate_dehydrogenase
from .glycolate_oxidase import add_glycolate_oxidase, add_glycolate_oxidase_yokota
from .glycolyl_coa_synthetase import (
    add_glycolyl_coa_synthetase,
    add_glycolyl_coa_synthetase_irrev,
)
from .glyoxylate_carboligase import add_glyoxylate_carboligase
from .glyoxylate_oxidase import add_glyoxylate_oxidase
from .hydroxyaspartate_hydrolase import add_hydroxyaspartate_hydrolase
from .hydroxypyruvate_isomerase import add_hydroxypyruvate_isomerase
from .lhc_deprotonation import add_lhc_deprotonation
from .lhc_protonation import add_lhc_protonation
from .lhc_state_transition import (
    add_state_transition_12,
    add_state_transition_21,
    add_state_transitions,
)
from .lumped_pgk_gadph import lumped_pgk_gadph
from .malate_dehydrogenase import add_malate_dehydrogenase
from .malate_synthase import add_malate_synthase
from .malic_enzyme import add_malic_enzyme
from .malony_coa_reductase import add_malonyl_coa_reductase
from .mda_reductase1 import add_mda_reductase1
from .mda_reductase2 import add_mda_reductase2
from .methylene_thf_dehydrogenase import add_methylene_thf_dehydrogenase
from .mthfc import add_mthfc
from .ndh import add_ndh
from .nitrogen_fixation import add_nitrogen_metabolism
from .oxalate_oxidase import add_oxalate_oxidase
from .oxaloacetate_formation import add_oxaloacetate_formation
from .oxidative_phosphorylation import add_oxidative_phosphorylation
from .pep_carboxylase import add_pep_carboxylase
from .phosphoglucomutase import add_phosphoglucomutase
from .phosphoglycerate_kinase import (
    add_phosphoglycerate_kinase,
    add_phosphoglycerate_kinase_poolman,
)
from .phosphoglycerate_mutase import add_phosphoglycerate_mutase
from .phosphoglycolate_phosphatase import (
    add_phosphoglycolate_influx,
    add_phosphoglycolate_phosphatase,
)
from .phosphoribulokinase import add_phosphoribulokinase
from .proton_leak import add_proton_leak
from .psi_psii import (
    add_energy_production,
    add_mehler,
    add_photosystems,
    add_ps2_cross_section,
    add_psi_2019,
    add_psi_2021,
    add_psii,
)
from .ptox import add_ptox
from .pyruvate_dehydrogenase import add_pyruvate_dehydrogenase
from .pyruvate_phosphate_dikinase import add_pyruvate_phosphate_dikinase
from .quencher import add_quenching_reaction
from .r1p_kinase import add_r1p_kinase
from .ribulose_phosphate_epimerase import add_ribulose_5_phosphate_3_epimerase
from .rpi import add_ribose_5_phosphate_isomerase
from .rubisco import add_rubisco, add_rubisco_poolman
from .sbpase import add_sbpase
from .serine_glyoxylate_transaminase import (
    add_serine_glyoxylate_transaminase,
    add_serine_glyoxylate_transaminase_irreversible,
)
from .succinyl_coa_synthetase import add_succinyl_coa_synthetase
from .tartronic_semialdehyde_reductase import add_tartronate_semialdehyde_reductase
from .tartronyl_coa_reductase import add_tartronyl_coa_reductase
from .thioesterase import add_thioesterase
from .thioredoxin import (
    add_cbb_pfd_linear_speedup,
    add_cbb_pfd_mm_speedup,
    add_e_relaxation,
    add_e_relaxation_2021,
    add_fd_tr_reductase,
    add_fd_tr_reductase_2021,
    add_nadph_tr_reductase,
    add_thioredoxin_regulation,
    add_thioredoxin_regulation2021,
    add_thioredoxin_regulation_from_nadph,
    add_tr_e_activation,
    add_tr_e_activation2021,
)
from .tp_export import (
    add_dhap_exporter,
    add_gap_exporter,
    add_pga_exporter,
    add_triose_phosphate_exporters,
)
from .transaldolase_f6p_gad_xyl_gap import add_transaldolase_f6p_gad_xyl_gap
from .transketolase_gad_s7p_eru_r5p import add_transketolase_gad_s7p_eru_r5p
from .transketolase_x5p_e4p_f6p_gap import add_transketolase_x5p_e4p_f6p_gap
from .transketolase_x5p_r5p_s7p_gap import add_transketolase_x5p_r5p_s7p_gap
from .triose_phosphate_isomerase import add_triose_phosphate_isomerase
from .violaxanthin_deepoxidase import add_violaxanthin_epoxidase
from .xylulose_kinase import add_xylulose_kinase
from .zeaxanthin_epoxidase import add_zeaxanthin_epoxidase

__all__ = [
    "add_acetyl_coa_acetyltransfer",
    "add_glycolaldehyde_dehydrogenase",
    "add_methylene_thf_dehydrogenase",
    "add_glucose_6_phosphate_isomerase_re",
    "lumped_pgk_gadph",
    "add_atp_consumption",
    "add_nadph_consumption",
    "add_hydroxyaspartate_aldolase",
    "add_thioesterase",
    "add_oxalate_oxidase",
    "add_ferredoxin_reductase",
    "add_ascorbate_peroxidase",
    "add_lhc_deprotonation",
    "add_mda_reductase2",
    "add_violaxanthin_epoxidase",
    "add_proton_leak",
    "add_sbpase",
    "add_quenching_reaction",
    "add_r1p_kinase",
    "add_hydroxyaspartate_hydrolase",
    "add_serine_glyoxylate_transaminase_irreversible",
    "add_serine_glyoxylate_transaminase",
    "add_glyoxylate_oxidase",
    "add_glutathion_reductase_irrev",
    "add_ribose_5_phosphate_isomerase",
    "add_malonyl_coa_reductase",
    "add_oxaloacetate_formation",
    "add_transaldolase_f6p_gad_xyl_gap",
    "add_erythrulose_kinase",
    "add_phosphoglycerate_kinase_poolman",
    "add_phosphoglycerate_kinase",
    "add_aspartate_aminotransferase",
    "add_ptox",
    "add_malic_enzyme",
    "add_aspartate_oxidoreductase",
    "add_atp_synthase_mmol_chl",
    "add_atp_synthase_mm",
    "add_atp_synthase_static_protons",
    "add_atp_synthase_energy_dependent",
    "add_atp_synthase_2024",
    "add_succinyl_coa_synthetase",
    "add_aspartate_nadp_oxidoreductase",
    "add_catalase",
    "add_glycine_decarboxylase_yokota",
    "add_glycine_decarboxylase_irreversible",
    "add_glycine_decarboxylase",
    "add_e4p_isomerase",
    "add_aldolase_dhap_gap_req",
    "add_formate_thf_ligase",
    "add_zeaxanthin_epoxidase",
    "add_malate_synthase",
    "add_tartronate_semialdehyde_reductase",
    "add_glycerate_kinase",
    "add_phosphoglycerate_mutase",
    "add_phosphoribulokinase",
    "add_mthfc",
    "add_ps2_cross_section",
    "add_psii",
    "add_psi_2019",
    "add_psi_2021",
    "add_mehler",
    "add_photosystems",
    "add_energy_production",
    "add_glycine_transaminase_yokota",
    "add_glycine_transaminase_irreversible",
    "add_glycine_transaminase",
    "add_ndh",
    "add_co2_dissolving",
    "add_aldolase_dhap_a5p_req",
    "add_mda_reductase1",
    "add_cyclic_electron_flow",
    "add_oxidative_phosphorylation",
    "add_hpa_outflux",
    "add_glycerate_dehydrogenase",
    "add_formate_dehydrogenase",
    "add_glycolate_dehydrogenase",
    "add_ribulose_5_phosphate_3_epimerase",
    "add_r1p_aldolase",
    "add_hydroxypyruvate_isomerase",
    "add_transketolase_gad_s7p_eru_r5p",
    "add_glutamate_dehydrogenase",
    "add_glycolate_oxidase_yokota",
    "add_glycolate_oxidase",
    "add_tartronyl_coa_reductase",
    "add_phosphoglycolate_influx",
    "add_phosphoglycolate_phosphatase",
    "add_pyruvate_phosphate_dikinase",
    "add_dehydroascorbate_reductase",
    "add_enolase",
    "add_glyoxylate_carboligase",
    "add_acetyl_coa_carboxyltransfer",
    "add_acetyl_coa_carboxyltransfer_1i",
    "add_e4p_epimerase",
    "add_glycolyl_coa_synthetase_irrev",
    "add_glycolyl_coa_synthetase",
    "add_cbb_pfd_linear_speedup",
    "add_cbb_pfd_mm_speedup",
    "add_fd_tr_reductase_2021",
    "add_fd_tr_reductase",
    "add_nadph_tr_reductase",
    "add_tr_e_activation",
    "add_tr_e_activation2021",
    "add_e_relaxation",
    "add_e_relaxation_2021",
    "add_thioredoxin_regulation",
    "add_thioredoxin_regulation2021",
    "add_thioredoxin_regulation_from_nadph",
    "add_a5p_aldolase",
    "add_pyruvate_dehydrogenase",
    "add_g1p_efflux",
    "add_malate_dehydrogenase",
    "add_phosphoglucomutase",
    "add_pga_exporter",
    "add_gap_exporter",
    "add_dhap_exporter",
    "add_triose_phosphate_exporters",
    "add_transketolase_x5p_e4p_f6p_gap",
    "add_a5p_isomerase",
    "add_pep_carboxylase",
    "add_triose_phosphate_isomerase",
    "add_gadph",
    "add_state_transition_12",
    "add_state_transition_21",
    "add_state_transitions",
    "add_lhc_protonation",
    "add_aldolase_dhap_e4p_req",
    "add_fnr_mmol_chl",
    "add_fnr_mm",
    "add_fnr_static",
    "add_fnr_energy_dependent",
    "add_nitrogen_metabolism",
    "add_carbonic_anhydrase_mass_action",
    "add_xylulose_kinase",
    "add_b6f",
    "add_b6f_2024",
    "add_fbpase",
    "add_rubisco_poolman",
    "add_rubisco",
    "add_acetoacetate_coa_ligase",
    "add_transketolase_x5p_r5p_s7p_gap",
]

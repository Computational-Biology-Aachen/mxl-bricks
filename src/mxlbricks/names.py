"""Standardized name constructors for metabolites, enzymes, and kinetic parameters."""

from __future__ import annotations

from typing import Literal

EMPTY: Literal[""] = ""

###############################################################################
# Parameter fns
###############################################################################


def dummy(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return standardized name for a dummy metabolite."""
    return loc("dummy", compartment, tissue)


def substrate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return standardized name for a generic substrate."""
    return loc("substrate", compartment, tissue)


def product(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return standardized name for a generic product."""
    return loc("product", compartment, tissue)


###############################################################################
# Parameter fns
###############################################################################


def loc(name: str, compartment: str = "", tissue: str = "") -> str:
    """Localise a component to a compartment and tissue."""
    return f"{name}{compartment}{tissue}"


def e0(enzyme: str) -> str:
    """Return standardized name for total enzyme concentration parameter."""
    return f"E0_{enzyme}"


def e(enzyme: str) -> str:
    """Return standardized name for enzyme state variable."""
    return f"E_{enzyme}"


def kcat(enzyme: str) -> str:
    """Return standardized name for catalytic rate constant."""
    return f"kcat_{enzyme}"


def vmax(enzyme: str) -> str:
    """Return standardized name for maximum reaction rate."""
    return f"vmax_{enzyme}"


def keq(enzyme: str) -> str:
    """Return standardized name for equilibrium constant."""
    return f"keq_{enzyme}"


def kre(enzyme: str) -> str:
    """Return standardized name for rapid-equilibrium rate constant."""
    return f"kre_{enzyme}"


def kf(enzyme: str) -> str:
    """Return standardized name for forward rate constant."""
    return f"kf_{enzyme}"


def kh(enzyme: str) -> str:
    """Hill constant."""
    return f"kh_{enzyme}"


def ksat(enzyme: str) -> str:
    """Return standardized name for saturation constant."""
    return f"ksat_{enzyme}"


def km(enzyme: str, substrate: str | None = None) -> str:
    """Return standardized name for Michaelis constant."""
    if substrate is None:
        return f"km_{enzyme}"
    return f"km_{enzyme}_{substrate}"


def kms(enzyme: str) -> str:
    """Return standardized name for Michaelis constant of substrate."""
    return km(enzyme, "s")


def kmp(enzyme: str) -> str:
    """Return standardized name for Michaelis constant of product."""
    return km(enzyme, "p")


def ki(enzyme: str, substrate: str | None = None) -> str:
    """Return standardized name for inhibition constant."""
    if substrate is None:
        return f"ki_{enzyme}"
    return f"ki_{enzyme}_{substrate}"


def ka(enzyme: str, substrate: str | None = None) -> str:
    """Return standardized name for activation constant."""
    if substrate is None:
        return f"ki_{enzyme}"
    return f"ki_{enzyme}_{substrate}"


###############################################################################
# Parameters / Variables
###############################################################################


def rt(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return dummy name for energy state."""
    return loc("RT", compartment, tissue)


def energy(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return dummy name for energy state."""
    return loc("energy", compartment, tissue)


def a0(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Photosystem II reaction center 0."""
    return loc("A0", compartment, tissue)


def a1(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Photosystem II reaction center 1."""
    return loc("A1", compartment, tissue)


def a2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Photosystem II reaction center 2."""
    return loc("A2", compartment, tissue)


def b0(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII reaction center state B0."""
    return loc("B0", compartment, tissue)


def b1(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII reaction center state B1."""
    return loc("B1", compartment, tissue)


def b2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII reaction center state B2."""
    return loc("B2", compartment, tissue)


def b3(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII reaction center state B3."""
    return loc("B3", compartment, tissue)


def ps2cs(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII cross-section."""
    return loc("PSII_cross_section", compartment, tissue)


def atp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Adenosine triphosphate."""
    return loc("ATP", compartment, tissue)


def adp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Adenosine diphosphate."""
    return loc("ADP", compartment, tissue)


def amp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Adenosine monophosphate."""
    return loc("AMP", compartment, tissue)


def nadph(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Nicotinamide adenine dinucleotide phosphate (reduced)."""
    return loc("NADPH", compartment, tissue)


def nadp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Nicotinamide adenine dinucleotide phosphate (oxidised)."""
    return loc("NADP", compartment, tissue)


def nadh(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Nicotinamide adenine dinucleotide (reduced)."""
    return loc("NADH", compartment, tissue)


def nad(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Nicotinamide adenine dinucleotide (oxidised)."""
    return loc("NAD", compartment, tissue)


def h(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Proton concentration."""
    return loc("protons", compartment, tissue)


def ph(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PH."""
    return loc("pH", compartment, tissue)


def pq_ox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Plastoquinone (oxidised)."""
    return loc("Plastoquinone (oxidised)", compartment, tissue)


def pq_red(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Plastoquinone (reduced)."""
    return loc("Plastoquinone (reduced)", compartment, tissue)


def pc_ox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Plastocyanine (oxidised)."""
    return loc("Plastocyanine (oxidised)", compartment, tissue)


def pc_red(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Plastocyanine (reduced)."""
    return loc("Plastocyanine (reduced)", compartment, tissue)


def fd_ox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ferredoxine (oxidised)."""
    return loc("Ferredoxine (oxidised)", compartment, tissue)


def fd_red(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ferredoxine (reduced)."""
    return loc("Ferredoxine (reduced)", compartment, tissue)


def lhc(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Light-harvesting complex."""
    return loc("Light-harvesting complex", compartment, tissue)


def lhcp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Protonated light-harvesting complex."""
    return loc("Light-harvesting complex (protonated)", compartment, tissue)


def psbs_de(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Deprotonated Psbs."""
    return loc("PsbS (de-protonated)", compartment, tissue)


def psbs_pr(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Protonated Psbs."""
    return loc("PsbS (protonated)", compartment, tissue)


def vx(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Violaxanthin."""
    return loc("Violaxanthin", compartment, tissue)


def zx(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Zeaxanthin."""
    return loc("Zeaxanthin", compartment, tissue)


def tr_ox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Thioredoxin (oxidised)."""
    return loc("Thioredoxin (oxidised)", compartment, tissue)


def tr_red(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Thioredoxin (reduced)."""
    return loc("Thioredoxin (reduced)", compartment, tissue)


def pi(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Orthophosphate."""
    return loc("Orthophosphate", compartment, tissue)


def pi_ext(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Orthophosphate (external)."""
    return loc("Orthophosphate (external)", compartment, tissue)


def ppi(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Diphosphate."""
    return loc("Diphosphate", compartment, tissue)


def co2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """CO2 (dissolved)."""
    return loc("CO2 (dissolved)", compartment, tissue)


def co2_atmosphere(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """CO2 (atmosphere)."""
    return loc("CO2 (atmosphere)", compartment, tissue)


def hco3(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """HCO3."""
    return loc("HCO3", compartment, tissue)


def o2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """O2 (dissolved)."""
    return loc("O2 (dissolved)", compartment, tissue)


def o2_atmosphere(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """O2 (atmosphere)."""
    return loc("O2 (atmosphere)", compartment, tissue)


def h2o2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """H2O2."""
    return loc("H2O2", compartment, tissue)


def pga(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """3PGA."""
    return loc("3PGA", compartment, tissue)


def pga2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """2PGA."""
    return loc("2PGA", compartment, tissue)


def pgo(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PGO."""
    return loc("PGO", compartment, tissue)


def bpga(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """BPGA."""
    return loc("BPGA", compartment, tissue)


def gap(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """GAP."""
    return loc("GAP", compartment, tissue)


def dhap(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """DHAP."""
    return loc("DHAP", compartment, tissue)


def fbp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """FBP."""
    return loc("FBP", compartment, tissue)


def f6p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """F6P."""
    return loc("F6P", compartment, tissue)


def g6p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """G6P."""
    return loc("G6P", compartment, tissue)


def g1p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """G1P."""
    return loc("G1P", compartment, tissue)


def sbp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """SBP."""
    return loc("SBP", compartment, tissue)


def s7p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """S7P."""
    return loc("S7P", compartment, tissue)


def erythrose(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Erythrose."""
    return loc("erythrose", compartment, tissue)


def e4p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Erythrose-4-phosphate."""
    return loc("E4P", compartment, tissue)


def erythrulose(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Erythrulose."""
    return loc("erythrulose", compartment, tissue)


def erythrulose_1p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """erythrulose_1p."""
    return loc("erythrulose_1p", compartment, tissue)


def erythrulose_4p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """erythrulose_4p."""
    return loc("erythrulose_4p", compartment, tissue)


def xylulose(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Xylulose."""
    return loc("xylulose", compartment, tissue)


def x5p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Xylulose-5-phosphate."""
    return loc("X5P", compartment, tissue)


def ribose(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ribose."""
    return loc("ribose", compartment, tissue)


def r1p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ribose-1-phosphate."""
    return loc("R1P", compartment, tissue)


def r5p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ribose-5-phosphate."""
    return loc("R5P", compartment, tissue)


def rubp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ribulose-1,5-bisphosphate."""
    return loc("RUBP", compartment, tissue)


def ru5p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """RU5P."""
    return loc("RU5P", compartment, tissue)


def o8p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """D-glycero-D-altro-octulose 8-phosphate."""
    # octulose 8-phosphate
    # D-glycero-D-altro-octulose 8-phosphate
    return loc("O8P", compartment, tissue)


def obp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """D-erythro-D-gluco-octose alpha-1,8-bisphosphate."""
    # octose 1,8-bisphosphate
    # D-erythro-D-gluco-octose α-1,8-bisphosphate  # noqa: RUF003
    return loc("OBP", compartment, tissue)


def starch(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Starch."""
    return loc("starch", compartment, tissue)


def mda(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """MDA."""
    return loc("MDA", compartment, tissue)


def dha(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """DHA."""
    return loc("DHA", compartment, tissue)


def ascorbate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ascorbate."""
    return loc("ascorbate", compartment, tissue)


def glutathion_red(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glutathion (reduced) / GSH."""
    return loc("GSH", compartment, tissue)


def glutathion_ox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glutathion (oxidised) / GSSG."""
    return loc("GSSG", compartment, tissue)


def glycine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glycine."""
    return loc("glycine", compartment, tissue)


def glyoxylate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glyoxylate."""
    return loc("glyoxylate", compartment, tissue)


def glycolate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glycolate."""
    return loc("glycolate", compartment, tissue)


def glycerate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glycerate."""
    return loc("glycerate", compartment, tissue)


def pyruvate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Pyruvate."""
    return loc("pyruvate", compartment, tissue)


def hydroxypyruvate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Hydroxypyruvate."""
    return loc("hydroxypyruvate", compartment, tissue)


def serine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Serine."""
    return loc("serine", compartment, tissue)


def pfd(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Photosynthetic Photon Flux Density."""
    return loc("PPFD", compartment, tissue)


def quencher(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Q."""
    return loc("Q", compartment, tissue)


def fluorescence(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Fluo."""
    return loc("Fluo", compartment, tissue)


def e_active(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """E_active."""
    return loc("E_active", compartment, tissue)


def e_inactive(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """E_inactive."""
    return loc("E_inactive", compartment, tissue)


def e_total(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """E_total."""
    return loc("E_total", compartment, tissue)


def oxoglutarate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """2-oxoglutarate."""
    return loc("2-oxoglutarate", compartment, tissue)


def glutamate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """GLU."""
    return loc("GLU", compartment, tissue)


def glutamine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """GLN."""
    return loc("GLN", compartment, tissue)


def nh4(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """NH4."""
    return loc("NH4", compartment, tissue)


def acetoacetate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """aka: ketobutyrate."""
    return loc("acetoacetate", compartment, tissue)


def coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """CoA."""
    return loc("CoA", compartment, tissue)


def acetyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """acetyl-CoA."""
    return loc("acetyl-CoA", compartment, tissue)


def acetoacetyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """acetoacetyl-CoA."""
    return loc("acetoacetyl-CoA", compartment, tissue)


def malonyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malonyl-CoA."""
    return loc("malonyl-CoA", compartment, tissue)


def formyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """formyl-CoA."""
    return loc("formyl-CoA", compartment, tissue)


def tartronyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """tartronyl-CoA."""
    return loc("tartronyl-CoA", compartment, tissue)


def succinyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """succinyl-CoA."""
    return loc("succinyl-CoA", compartment, tissue)


def glycolyl_coa(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolyl-CoA."""
    return loc("glycolyl-CoA", compartment, tissue)


def malonate_s_aldehyde(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malonate-s-aldehyde."""
    return loc("malonate-s-aldehyde", compartment, tissue)


def succinate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Succinate."""
    return loc("succinate", compartment, tissue)


def formate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Return formate name."""
    return loc("formate", compartment, tissue)


def thf(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """THF."""
    return loc("THF", compartment, tissue)


def formyl_thf(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """10-formyl-THF."""
    return loc("10-formyl-THF", compartment, tissue)


def methylene_thf(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """methylene-THF."""
    return loc("methylene-THF", compartment, tissue)


def methenyl_thf(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """methenyl-THF."""
    return loc("methenyl-THF", compartment, tissue)


def aspartate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Aspartate."""
    return loc("aspartate", compartment, tissue)


def hydroxyaspartate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Hydroxyaspartate."""
    return loc("hydroxyaspartate", compartment, tissue)


def iminoaspartate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Iminoaspartate."""
    return loc("iminoaspartate", compartment, tissue)


def oxalate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Oxalate (oxalic acid / ethanedioic acid)."""
    # synonyms: ethanedoic acid, oxalic acid
    return loc("oxalate", compartment, tissue)


def oxaloacetate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Oxaloacetate."""
    return loc("oxaloacetate", compartment, tissue)


def malate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Malate."""
    return loc("malate", compartment, tissue)


def pep(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Phosphoenolpyruvat."""
    return loc("PEP", compartment, tissue)


def tartronate_semialdehyde(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """tartronate_semialdehyde."""
    return loc("tartronate_semialdehyde", compartment, tissue)


def arabinose_5_phosphate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """A5P."""
    return loc("A5P", compartment, tissue)


def glycolaldehyde(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glycolaldehyde."""
    return loc("glycolaldehyde", compartment, tissue)


def alanine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Alanine."""
    return loc("alanine", compartment, tissue)


def arginine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Arginine."""
    return loc("arginine", compartment, tissue)


def asparagine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Asparagine."""
    return loc("asparagine", compartment, tissue)


def cysteine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Cysteine."""
    return loc("cysteine", compartment, tissue)


def histidine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Histidine."""
    return loc("histidine", compartment, tissue)


def isoleucine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Isoleucine."""
    return loc("isoleucine", compartment, tissue)


def leucine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Leucine."""
    return loc("leucine", compartment, tissue)


def lysine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Lysine."""
    return loc("lysine", compartment, tissue)


def methionine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Methionine."""
    return loc("methionine", compartment, tissue)


def phenylalanine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Phenylalanine."""
    return loc("phenylalanine", compartment, tissue)


def proline(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Proline."""
    return loc("proline", compartment, tissue)


def threonine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Threonine."""
    return loc("threonine", compartment, tissue)


def tryptophan(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Tryptophan."""
    return loc("tryptophan", compartment, tissue)


def tyrosine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Tyrosine."""
    return loc("tyrosine", compartment, tissue)


def valine(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Valine."""
    return loc("valine", compartment, tissue)


###############################################################################
# Moieties
###############################################################################


def total_adenosines(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """A*P."""
    return loc("A*P", compartment, tissue)


def total_nadp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """NADP*."""
    return loc("NADP*", compartment, tissue)


def total_nad(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """NAD*."""
    return loc("NAD*", compartment, tissue)


def total_pq(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PQ_tot."""
    return loc("PQ_tot", compartment, tissue)


def total_pc(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PC_tot."""
    return loc("PC_tot", compartment, tissue)


def total_ascorbate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ASC_tot*."""
    return loc("ASC_tot*", compartment, tissue)


def total_ferredoxin(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Fd*."""
    return loc("Fd*", compartment, tissue)


def total_glutamate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glu+Oxo."""
    return loc("Glu+Oxo", compartment, tissue)


def total_glutathion(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Glutathion_tot."""
    return loc("Glutathion_tot", compartment, tissue)


def total_lhc(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """LHC_tot."""
    return loc("LHC_tot", compartment, tissue)


def total_orthophosphate(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Pi_tot."""
    return loc("Pi_tot", compartment, tissue)


def total_carotenoids(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Carotenoids_tot."""
    return loc("Carotenoids_tot", compartment, tissue)


def total_thioredoxin(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Thioredoxin_tot."""
    return loc("Thioredoxin_tot", compartment, tissue)


def total_psbs(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSBS_tot."""
    return loc("PSBS_tot", compartment, tissue)


###############################################################################
# Reactions / Enzymes
###############################################################################


def aspartate_aminotransferase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """aspartate_aminotransferase."""
    return loc("aspartate_aminotransferase", compartment, tissue)


def aspartate_oxidoreductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """aspartate_oxidoreductase."""
    return loc("aspartate_oxidoreductase", compartment, tissue)


def oxidative_phosphorylation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """oxidative_phosphorylation."""
    return loc("oxidative_phosphorylation", compartment, tissue)


def oxalate_oxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """oxalate_oxidase."""
    return loc("oxalate_oxidase", compartment, tissue)


def hydroxypyruvate_isomerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """hydroxypyruvate_isomerase."""
    return loc("hydroxypyruvate_isomerase", compartment, tissue)


def pgk_gadph(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """pgk_gadph."""
    return loc("pgk_gadph", compartment, tissue)


def glycolaldehyde_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolaldehyde_dehydrogenase."""
    return loc("glycolaldehyde_dehydrogenase", compartment, tissue)


def a5p_aldolase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """a5p_aldolase."""
    return loc("a5p_aldolase", compartment, tissue)


def a5p_isomerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """a5p_isomerase."""
    return loc("a5p_isomerase", compartment, tissue)


def r1p_aldolase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """r1p_aldolase."""
    return loc("r1p_aldolase", compartment, tissue)


def r1p_kinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """r1p_kinase."""
    return loc("r1p_kinase", compartment, tissue)


def transaldolase_f6p_gad_gap_xyl(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """transaldolase_f6p_gad_gap_xyl."""
    return loc("transaldolase_f6p_gad_gap_xyl", compartment, tissue)


def transketolase_gad_s7p_r5p_eru(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """transketolase_gad_s7p_r5p_eru."""
    return loc("transketolase_gad_s7p_r5p_eru", compartment, tissue)


def transketolase_f6p_gad_gap_xylulose(
    compartment: str = EMPTY, tissue: str = EMPTY
) -> str:
    """transketolase_f6p_gad_gap_xylulose."""
    return loc("transketolase_f6p_gad_gap_xylulose", compartment, tissue)


def enolase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Enolase."""
    return loc("enolase", compartment, tissue)


def erythrulose_kinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """erythrulose_kinase."""
    return loc("erythrulose_kinase", compartment, tissue)


def e4p_epimerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """e4p_epimerase."""
    return loc("e4p_epimerase", compartment, tissue)


def e4p_isomerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """e4p_isomerase."""
    return loc("e4p_isomerase", compartment, tissue)


def xylulose_kinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """xylulose_kinase."""
    return loc("xylulose_kinase", compartment, tissue)


def acetoacetate_coa_ligase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """acetoacetate_coa_ligase."""
    return loc("acetoacetate_coa_ligase", compartment, tissue)


def acetyl_coa_acetyltransfer(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """acetyl_coa_acetyltransfer."""
    return loc("acetyl_coa_acetyltransfer", compartment, tissue)


def acetyl_coa_carboxyltransfer(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """acetyl_coa_carboxyltransfer."""
    return loc("acetyl_coa_carboxyltransfer", compartment, tissue)


def aldolase_dhap_gap(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """aldolase_dhap_gap."""
    return loc("aldolase_dhap_gap", compartment, tissue)


def aldolase_dhap_e4p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """aldolase_dhap_e4p."""
    return loc("aldolase_dhap_e4p", compartment, tissue)


def ascorbate_peroxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ascorbate_peroxidase."""
    return loc("ascorbate_peroxidase", compartment, tissue)


def atp_synthase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """atp_synthase."""
    return loc("atp_synthase", compartment, tissue)


def b6f(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """b6f."""
    return loc("b6f", compartment, tissue)


def bkace(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Bkace."""
    return loc("bkace", compartment, tissue)


def carbonic_anhydrase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """carbonic_anhydrase."""
    return loc("carbonic_anhydrase", compartment, tissue)


def catalase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Catalase."""
    return loc("catalase", compartment, tissue)


def cyclic_electron_flow(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """cyclic_electron_flow."""
    return loc("cyclic_electron_flow", compartment, tissue)


def co2_dissolving(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """co2_dissolving."""
    return loc("co2_dissolving", compartment, tissue)


def coa_transf_a(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """coa_transf_a."""
    return loc("coa_transf_a", compartment, tissue)


def coa_transf_b(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """coa_transf_b."""
    return loc("coa_transf_b", compartment, tissue)


def dehydroascorbate_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """dehydroascorbate_reductase."""
    return loc("dehydroascorbate_reductase", compartment, tissue)


def ex_atp(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_atp."""
    return loc("ex_atp", compartment, tissue)


def ex_nadph(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_nadph."""
    return loc("ex_nadph", compartment, tissue)


def ex_g1p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_g1p."""
    return loc("ex_g1p", compartment, tissue)


def ex_pga(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_pga."""
    return loc("ex_pga", compartment, tissue)


def ex_gap(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_gap."""
    return loc("ex_gap", compartment, tissue)


def ex_dhap(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ex_dhap."""
    return loc("ex_dhap", compartment, tissue)


def fbpase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Fbpase."""
    return loc("fbpase", compartment, tissue)


def ferredoxin_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ferredoxin_reductase."""
    return loc("ferredoxin_reductase", compartment, tissue)


def ferredoxin_thioredoxin_reductase(
    compartment: str = EMPTY, tissue: str = EMPTY
) -> str:
    """ferredoxin_thioredoxin_reductase."""
    return loc("ferredoxin_thioredoxin_reductase", compartment, tissue)


def nadph_thioredoxin_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """nadph_thioredoxin_reductase."""
    return loc("nadph_thioredoxin_reductase", compartment, tissue)


def fnr(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Fnr."""
    return loc("fnr", compartment, tissue)


def formate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """formate_dehydrogenase."""
    return loc("formate_dehydrogenase", compartment, tissue)


def formate_thf_ligase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """formate_thf_ligase."""
    return loc("formate_thf_ligase", compartment, tissue)


def g6pi(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """g6pi."""
    return loc("g6pi", compartment, tissue)


def gadph(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Gadph."""
    return loc("gadph", compartment, tissue)


def glutathion_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glutathion_reductase."""
    return loc("glutathion_reductase", compartment, tissue)


def glycerate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycerate_dehydrogenase."""
    return loc("glycerate_dehydrogenase", compartment, tissue)


def glycerate_kinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycerate_kinase."""
    return loc("glycerate_kinase", compartment, tissue)


def glycine_decarboxylase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycine_decarboxylase."""
    return loc("glycine_decarboxylase", compartment, tissue)


def glycine_hydroxymethyltransferase(
    compartment: str = EMPTY, tissue: str = EMPTY
) -> str:
    """glycine_hydroxymethyltransferase."""
    return loc("glycine_hydroxymethyltransferase", compartment, tissue)


def glycine_transaminase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycine_transaminase."""
    return loc("glycine_transaminase", compartment, tissue)


def glycolate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolate_dehydrogenase."""
    return loc("glycolate_dehydrogenase", compartment, tissue)


def glycolate_oxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolate_oxidase."""
    return loc("glycolate_oxidase", compartment, tissue)


def glyoxylate_oxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glyoxylate_oxidase."""
    return loc("glyoxylate_oxidase", compartment, tissue)


def glycolyl_coa_carboxylase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolyl_coa_carboxylase."""
    return loc("glycolyl_coa_carboxylase", compartment, tissue)


def glycolyl_coa_synthetase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glycolyl_coa_synthetase."""
    return loc("glycolyl_coa_synthetase", compartment, tissue)


def glyoxylate_carboligase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glyoxylate_carboligase."""
    return loc("glyoxylate_carboligase", compartment, tissue)


def hydroxyaspartate_aldolase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """hydroxyaspartate_aldolase."""
    return loc("hydroxyaspartate_aldolase", compartment, tissue)


def hydroxyaspartate_hydrolase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """hydroxyaspartate_hydrolase."""
    return loc("hydroxyaspartate_hydrolase", compartment, tissue)


def lhc_deprotonation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """lhc_deprotonation."""
    return loc("lhc_deprotonation", compartment, tissue)


def lhc_protonation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """lhc_protonation."""
    return loc("lhc_protonation", compartment, tissue)


def lhc_state_transition_12(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """lhc_state_transition_12."""
    return loc("lhc_state_transition_12", compartment, tissue)


def lhc_state_transition_21(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """lhc_state_transition_21."""
    return loc("lhc_state_transition_21", compartment, tissue)


def malate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malate_dehydrogenase."""
    return loc("malate_dehydrogenase", compartment, tissue)


def malate_synthase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malate_synthase."""
    return loc("malate_synthase", compartment, tissue)


def malic_enzyme(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malic_enzyme."""
    return loc("malic_enzyme", compartment, tissue)


def malonyl_coa_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """malonyl_coa_reductase."""
    return loc("malonyl_coa_reductase", compartment, tissue)


def mda_reductase1(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Monodehydroascorbate reductase 1."""
    # FIXME
    return loc("mda_reductase_1", compartment, tissue)


def mda_reductase2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Monodehydroascorbate reductase 2."""
    # FIXME
    return loc("mda_reductase_2", compartment, tissue)


def mehler(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Mehler."""
    return loc("mehler", compartment, tissue)


def methylene_thf_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """methylene_thf_dehydrogenase."""
    return loc("methylene_thf_dehydrogenase", compartment, tissue)


def ndh(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Ndh."""
    return loc("ndh", compartment, tissue)


def mthfc(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Methenyltetrahydrofolate cyclohydrolase."""
    # FIXME: look up proper name of this
    return loc("mthfc", compartment, tissue)


def nitrogen_fixation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """nitrogen_fixation."""
    return loc("nitrogen_fixation", compartment, tissue)


def oxaloacetate_formation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """oxaloacetate_formation."""
    return loc("oxaloacetate_formation", compartment, tissue)


def pep_carboxylase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """pep_carboxylase."""
    return loc("pep_carboxylase", compartment, tissue)


def phosphoglucomutase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Phosphoglucomutase."""
    return loc("phosphoglucomutase", compartment, tissue)


def phosphoglycerate_kinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """phosphoglycerate_kinase."""
    return loc("phosphoglycerate_kinase", compartment, tissue)


def phosphoglycerate_mutase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """phosphoglycerate_mutase."""
    return loc("phosphoglycerate_mutase", compartment, tissue)


def phosphoglycolate_phosphatase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """phosphoglycolate_phosphatase."""
    return loc("phosphoglycolate_phosphatase", compartment, tissue)


def phosphoribulokinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Phosphoribulokinase."""
    return loc("phosphoribulokinase", compartment, tissue)


def proton_leak(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """proton_leak."""
    return loc("proton_leak", compartment, tissue)


def petc(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PETC."""
    return loc("PETC", compartment, tissue)


def ps1(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSI."""
    return loc("PSI", compartment, tissue)


def ps2(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PSII."""
    return loc("PSII", compartment, tissue)


def ptox(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """PTOX."""
    return loc("PTOX", compartment, tissue)


def pyruvate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """pyruvate_dehydrogenase."""
    return loc("pyruvate_dehydrogenase", compartment, tissue)


def pyruvate_phosphate_dikinase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """pyruvate_phosphate_dikinase."""
    return loc("pyruvate_phosphate_dikinase", compartment, tissue)


def ribulose_phosphate_epimerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ribulose_phosphate_epimerase."""
    return loc("ribulose_phosphate_epimerase", compartment, tissue)


def ribose_phosphate_isomerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """ribose_phosphate_isomerase."""
    return loc("ribose_phosphate_isomerase", compartment, tissue)


def rubisco(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Rubisco."""
    return loc("rubisco", compartment, tissue)


def rubisco_carboxylase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """rubisco_carboxylase."""
    return loc("rubisco_carboxylase", compartment, tissue)


def rubisco_oxygenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """rubisco_oxygenase."""
    return loc("rubisco_oxygenase", compartment, tissue)


def sbpase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """SBPase."""
    return loc("SBPase", compartment, tissue)


def serine_glyoxylate_transaminase(
    compartment: str = EMPTY, tissue: str = EMPTY
) -> str:
    """serine_glyoxylate_transaminase."""
    return loc("serine_glyoxylate_transaminase", compartment, tissue)


def succinyl_coa_synthetase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """succinyl_coa_synthetase."""
    return loc("succinyl_coa_synthetase", compartment, tissue)


def tartronate_semialdehyde_reductase(
    compartment: str = EMPTY, tissue: str = EMPTY
) -> str:
    """tartronate_semialdehyde_reductase."""
    return loc("tartronate_semialdehyde_reductase", compartment, tissue)


def tartronyl_coa_reductase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """tartronyl_coa_reductase."""
    return loc("tartronyl_coa_reductase", compartment, tissue)


def thioesterase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Thioesterase."""
    return loc("thioesterase", compartment, tissue)


def transketolase_gap_f6p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """transketolase_gap_f6p."""
    return loc("transketolase_gap_f6p", compartment, tissue)


def transketolase_gap_s7p(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """transketolase_gap_s7p."""
    return loc("transketolase_gap_s7p", compartment, tissue)


def triose_phosphate_isomerase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """triose_phosphate_isomerase."""
    return loc("triose_phosphate_isomerase", compartment, tissue)


def violaxanthin_deepoxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """violaxanthin_deepoxidase."""
    return loc("violaxanthin_deepoxidase", compartment, tissue)


def zeaxanthin_epoxidase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """zeaxanthin_epoxidase."""
    return loc("zeaxanthin_epoxidase", compartment, tissue)


def gogat(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """Gogat."""
    return loc("gogat", compartment, tissue)


def glutamine_synthase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glutamine_synthase."""
    return loc("glutamine_synthase", compartment, tissue)


def glutamate_dehydrogenase(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """glutamate_dehydrogenase."""
    return loc("glutamate_dehydrogenase", compartment, tissue)


def light_speedup(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """light_speedup."""
    return loc("light_speedup", compartment, tissue)


def tr_activation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """tr_activation."""
    return loc("tr_activation", compartment, tissue)


def tr_inactivation(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
    """tr_inactivation."""
    return loc("tr_inactivation", compartment, tissue)


def convf() -> str:
    """Conversion factor name."""
    return "convf"

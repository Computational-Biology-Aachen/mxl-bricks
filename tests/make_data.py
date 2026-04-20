# %%
import numpy as np
import pandas as pd

from mxlbricks import (
    get_matuszynska2016_npq,
    get_matuszynska2016_phd,
    get_matuszynska2019,
    get_poolman2000,
    get_saadat2021,
    get_yokota1985,
    names,
)

# %%
model = get_matuszynska2016_npq(mode="matrix")
y0 = {
    "pq_red": 0,
    "protons": 6.32975752e-05,
    "vmax_atp_synthase": 0,
    "atp": 25.0,
    "psbs_de": 1,
    "vx": 1,
}

pd.DataFrame(
    {
        ppfd: model.update_parameter(names.pfd(), ppfd).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for ppfd in np.arange(100, 1500, 100, dtype=float)
    }
).T.to_csv("assets/matuszynska-2016-npq-args-by-ppfd.csv")

# %%
model = get_matuszynska2016_phd(mode="matrix")
y0 = {
    "ATP": 1.6999999999999997,
    "Plastoquinone (oxidised)": 4.706348349506148,
    "Plastocyanine (oxidised)": 3.9414515288091567,
    "Ferredoxine (oxidised)": 3.7761613271207324,
    "protons_lumen": 7.737821100836988,
    "Light-harvesting complex": 0.5105293511676007,
    "PsbS (de-protonated)": 0.5000000001374878,
    "Violaxanthin": 0.09090909090907397,
}

pd.DataFrame(
    {
        ppfd: model.update_parameter(names.pfd(), ppfd).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for ppfd in np.arange(100, 1500, 100, dtype=float)
    }
).T.to_csv("assets/matuszynska-2016-phd-args-by-ppfd.csv")

# %%
model = get_matuszynska2019(mode="matrix")
y0 = {
    "3PGA": 0.9928653922138561,
    "BPGA": 0.0005297732935310749,
    "GAP": 0.0062663539939955834,
    "DHAP": 0.13785977143668732,
    "FBP": 0.006133532145409954,
    "F6P": 0.31271973359685457,
    "G6P": 0.719255387166192,
    "G1P": 0.041716812452951633,
    "SBP": 0.013123745088361893,
    "S7P": 0.15890073845176905,
    "E4P": 0.007322797350442026,
    "X5P": 0.022478763225333428,
    "R5P": 0.037651927659696716,
    "RUBP": 0.13184790283048484,
    "RU5P": 0.015060770937455408,
    "ATP": 1.612922506604933,
    "Ferredoxine (oxidised)": 3.8624032084329674,
    "protons_lumen": 0.002208423037307405,
    "Light-harvesting complex": 0.80137477470646,
    "NADPH": 0.491395685599137,
    "Plastocyanine (oxidised)": 1.885391998090184,
    "Plastoquinone (oxidised)": 10.991562708096392,
    "PsbS (de-protonated)": 0.9610220887579118,
    "Violaxanthin": 0.9514408605906095,
}

pd.DataFrame(
    {
        ppfd: model.update_parameter(names.pfd(), ppfd).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for ppfd in np.arange(100, 1500, 100, dtype=float)
    }
).T.to_csv("assets/matuszynska-2019-args-by-ppfd.csv")

# %%
model = get_saadat2021(mode="matrix")
y0 = {
    "3PGA": 0.9167729479368978,
    "BPGA": 0.0003814495319659031,
    "GAP": 0.00580821050261484,
    "DHAP": 0.1277806166216142,
    "FBP": 0.005269452472931973,
    "F6P": 0.2874944558066638,
    "G6P": 0.6612372482712676,
    "G1P": 0.03835176039761378,
    "SBP": 0.011101373736607443,
    "S7P": 0.1494578301900007,
    "E4P": 0.00668295494870102,
    "X5P": 0.020988553174809618,
    "R5P": 0.035155825913785584,
    "RUBP": 0.11293260727162346,
    "RU5P": 0.014062330254191594,
    "ATP": 1.4612747767895344,
    "Ferredoxine (oxidised)": 3.715702384326767,
    "protons_lumen": 0.002086128887296243,
    "Light-harvesting complex": 0.7805901436176024,
    "NADPH": 0.5578718406315588,
    "Plastocyanine (oxidised)": 1.8083642974980014,
    "Plastoquinone (oxidised)": 10.251099271612473,
    "PsbS (de-protonated)": 0.9667381262477079,
    "Violaxanthin": 0.9629870646993118,
    "MDA": 2.0353396709300447e-07,
    "H2O2": 1.2034405327140102e-07,
    "DHA": 1.0296456279861962e-11,
    "GSSG": 4.99986167652437e-12,
    "Thioredoxin (oxidised)": 0.9334426859846461,
    "E_inactive": 3.6023635680406634,
}

pd.DataFrame(
    {
        ppfd: model.update_parameter(names.pfd(), ppfd).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for ppfd in np.arange(100, 1500, 100, dtype=float)
    }
).T.to_csv("assets/saadat-2021-args-by-ppfd.csv")

# %%
model = get_poolman2000()
y0 = {
    "3PGA": 0.6387788347932627,
    "BPGA": 0.0013570885908749779,
    "GAP": 0.011259431827358068,
    "DHAP": 0.24770748227012374,
    "FBP": 0.01980222074817044,
    "F6P": 1.093666906864421,
    "G6P": 2.5154338857582377,
    "G1P": 0.14589516537322303,
    "SBP": 0.09132688566151095,
    "S7P": 0.23281380022778891,
    "E4P": 0.02836065066520614,
    "X5P": 0.03647242425941113,
    "R5P": 0.06109130988031577,
    "RUBP": 0.2672164362349537,
    "RU5P": 0.0244365238237522,
    "ATP": 0.43633201706180874,
}


pd.DataFrame(
    {
        i: model.update_parameter(names.e0(names.rubisco_carboxylase()), i).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for i in np.arange(0.1, 2.0, 0.1, dtype=float).round(1)
    }
).T.to_csv("assets/poolman-2000-args-by-influx.csv")

# %%
model = get_yokota1985()
y0 = {
    "glycolate": 0.09,
    "glyoxylate": 0.7964601770483386,
    "glycine": 8.999999999424611,
    "serine": 2.5385608670239126,
    "hydroxypyruvate": 0.009782608695111009,
    "H2O2": 0.010880542843616855,
}

pd.DataFrame(
    {
        i: model.update_parameter(
            names.kf(names.phosphoglycolate_phosphatase()), i
        ).get_args(
            y0,
            include_time=False,
            include_derived_parameters=True,
            include_derived_variables=True,
            include_parameters=False,
            include_reactions=True,
            include_readouts=True,
            include_surrogate_fluxes=True,
            include_surrogate_variables=True,
        )
        for i in np.arange(10, 100, 10, dtype=float)
    }
).T.to_csv("assets/yokota-1985-args-by-influx.csv")


# %%

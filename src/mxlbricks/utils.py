"""Utility helpers for building mxlpy models from mxlbricks."""

from collections.abc import Callable, Mapping

import numpy as np
import pandas as pd
import sympy
from mxlpy import Derived, Model
from mxlpy.fns import mul

from mxlbricks import names as n
from mxlbricks.fns import mass_action_1s


def static(
    model: Model,
    name: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Add a fixed parameter to model and return its name."""
    model.add_parameter(name, value, unit=unit, source=source)
    return name


def fcbb_regulated(model: Model, name: str, value: float) -> str:
    """Add a parameter scaled by the CBB light-speedup factor and return its derived name."""
    new_name = f"{name}_fcbb"

    model.add_parameter(name, value)
    model.add_derived(new_name, mul, args=[name, n.light_speedup()])
    return new_name


def thioredixon_regulated(model: Model, name: str, value: float) -> str:
    """Add a parameter scaled by the active enzyme fraction and return its derived name."""
    new_name = f"{name}_active"
    model.add_parameter(name, value)
    model.add_derived(new_name, mul, args=[name, n.e_active()])
    return new_name


def filter_stoichiometry(
    model: Model,
    stoichiometry: Mapping[str, float | Derived],
    optional: dict[str, float] | None = None,
) -> Mapping[str, float | Derived]:
    """Only use components that are actually compounds in the model."""
    variables = model.get_raw_variables(as_copy=False)

    new = {}
    for k, v in stoichiometry.items():
        if k in variables:
            new[k] = v
        elif k not in model._ids:  # noqa: SLF001
            msg = f"Missing component {k}"
            raise KeyError(msg)

    optional = {} if optional is None else optional
    new |= {k: v for k, v in optional.items() if k in variables}
    return new


def default_name(name: str | None, name_fn: Callable[[], str]) -> str:
    """Return name if given, otherwise call name_fn to get the default."""
    if name is None:
        return name_fn()
    return name


def default_par(
    model: Model,
    *,
    par: str | None,
    name: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return par if given, otherwise add a static parameter with the default value."""
    if par is not None:
        return par
    return static(model=model, name=name, value=value, unit=unit, source=source)


def default_keq(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the equilibrium constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.keq(rxn),
    )


def default_kf(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the forward rate constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.kf(rxn),
    )


def default_km(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    subs: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the Michaelis constant parameter for a substrate."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.km(rxn, subs),
    )


def default_kms(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the substrate Michaelis constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.kms(rxn),
    )


def default_kmp(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the product Michaelis constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.kmp(rxn),
    )


def default_ki(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the inhibition constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.ki(rxn),
    )


def default_kis(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    substrate: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the substrate-specific inhibition constant parameter."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.ki(rxn, substrate),
    )


def default_kre(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the rapid-equilibrium rate constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.kre(rxn),
    )


def default_e0(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the total enzyme concentration parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.e0(rxn),
    )


def default_kcat(
    model: Model,
    *,
    par: str | None,
    rxn: str,
    value: float,
    unit: sympy.Expr | None = None,
    source: str | None = None,
) -> str:
    """Return or create the catalytic rate constant parameter for a reaction."""
    return default_par(
        model=model,
        par=par,
        value=value,
        unit=unit,
        source=source,
        name=n.kcat(rxn),
    )


def default_vmax(
    model: Model,
    *,
    e0: str | None,
    kcat: str | None,
    rxn: str,
    e0_value: float,
    e0_unit: sympy.Expr | None = None,
    e0_source: str | None = None,
    kcat_value: float,
    kcat_unit: sympy.Expr | None = None,
    kcat_source: str | None = None,
) -> str:
    """Add vmax = kcat * E0 as a derived quantity and return its name."""
    e0 = default_e0(
        model=model,
        par=e0,
        rxn=rxn,
        value=e0_value,
        unit=e0_unit,
        source=e0_source,
    )
    kcat = default_kcat(
        model=model,
        par=kcat,
        rxn=rxn,
        value=kcat_value,
        unit=kcat_unit,
        source=kcat_source,
    )
    model.add_derived(vmax := n.vmax(rxn), fn=mass_action_1s, args=[kcat, e0])
    return vmax


def get_idxs_of_peaks(light: pd.Series) -> list[list[int]]:
    """Return groups of indices where light equals its maximum value."""
    cnt = 0
    idx_groups = []

    arr = light.to_numpy()
    max_light = light.max()
    while cnt < len(arr):
        if arr[cnt] == max_light:
            # temporary container for all F==maxlight. For each peak it is renewed
            idxs = []
            while cnt != len(arr) and arr[cnt] == max_light:
                idxs.append(cnt)
                cnt += 1
            idx_groups.append(idxs)
        else:
            cnt += 1
    return idx_groups


def get_idxs_of_peaks_(light: pd.Series) -> list[pd.Index]:
    """Return peak indices as pandas Index objects instead of integer lists."""
    mask = light == max(light)
    return list(mask[mask].groupby((1 - mask).cumsum()).apply(lambda x: x.index).values)


def get_npq(ppfd: pd.Series, fluorescence: pd.Series) -> pd.Series:
    """Calculate non-photochemical quenching from PAM simulation peak points.

    Returns
    -------
    Fm: Fm (first element of list) and Fm' values
    NPQ: Calculated NPQ values
    tm: Exact time points of peaks in PAM trace
    Fo: Fo (first element of list) and Ft' values
    to: Exact time points of Fo and Ft' values

    """
    # container for lists. Each list contains the positions of fluorescence values for one peak
    # container for position of Fo'
    t = ppfd.index.to_numpy()
    fluo = fluorescence.to_numpy()

    # Fm is the maximal value for each peak sequence
    max_fluo_per_group = [i[np.argmax(fluo[i])] for i in get_idxs_of_peaks(ppfd)]
    index = t[max_fluo_per_group]
    Fm = fluo[max_fluo_per_group]
    return pd.Series((Fm[0] - Fm) / Fm, name="NPQ", index=index)

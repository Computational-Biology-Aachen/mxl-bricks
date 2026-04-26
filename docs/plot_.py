"""Plotting helpers for mxlbricks model outputs."""

from __future__ import annotations

import pandas as pd
from mxlpy import plot
from mxlpy.plot import (
    Color,
    FigAxs,
    Linestyle,
    _partition_by_order_of_magnitude,
    _split_large_groups,
)


def _combine_small_groups(
    groups: list[list[str]], min_group_size: int
) -> list[list[str]]:
    """Combine adjacent groups if their size is smaller than min_group_size.

    Args:
        groups: List of lists of strings
        min_group_size: Minimum size for a group to remain separate

    Returns:
        List of lists with small adjacent groups combined

    """
    result = []
    current_group = groups[0]

    for next_group in groups[1:]:
        if len(current_group) < min_group_size:
            current_group.extend(next_group)
        else:
            result.append(current_group)
            current_group = next_group

    # Last group
    if len(current_group) < min_group_size:
        result[-1].extend(current_group)
    else:
        result.append(current_group)
    return result


def line_autogrouped(
    s: pd.Series | pd.DataFrame,
    *,
    n_cols: int = 2,
    col_width: float = 4,
    row_height: float = 3,
    min_group_size: int = 1,
    max_group_size: int = 6,
    grid: bool = True,
    xlabel: str | None = None,
    ylabel: str | None = None,
    color: Color | list[list[Color]] | None = None,
    linewidth: float | None = None,
    linestyle: Linestyle | None = None,
) -> FigAxs:
    """Plot a Series or DataFrame with auto-grouped subplots by order of magnitude."""
    group_names = (
        _partition_by_order_of_magnitude(s)
        if isinstance(s, pd.Series)
        else _partition_by_order_of_magnitude(s.max())
    )
    group_names = _combine_small_groups(group_names, min_group_size=min_group_size)
    group_names = _split_large_groups(group_names, max_size=max_group_size)

    groups: list[pd.Series] | list[pd.DataFrame] = (
        [s.loc[group] for group in group_names]
        if isinstance(s, pd.Series)
        else [s.loc[:, group] for group in group_names]
    )

    return plot.lines_grouped(
        groups,
        n_cols=n_cols,
        col_width=col_width,
        row_height=row_height,
        grid=grid,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        xlabel=xlabel,
        ylabel=ylabel,
    )

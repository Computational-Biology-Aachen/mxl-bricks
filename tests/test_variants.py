import pandas as pd

from mxlbricks import (
    get_matuszynska2016_npq,
    get_matuszynska2016_phd,
    get_matuszynska2019,
    get_saadat2021,
)
from mxlbricks import names as n


def test_variants_2016_npq() -> None:
    a1 = get_matuszynska2016_npq(mode="matrix").get_args()
    a2 = get_matuszynska2016_npq(mode="analytical").get_args()
    a3 = get_matuszynska2016_npq(mode="analytical-split").get_args()

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2016_phd() -> None:
    a1 = get_matuszynska2016_phd(mode="matrix").get_args()
    a2 = get_matuszynska2016_phd(mode="analytical").get_args()
    a3 = get_matuszynska2016_phd(mode="analytical-split").get_args()

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2019() -> None:
    a1 = get_matuszynska2019(mode="matrix").get_args()
    a2 = get_matuszynska2019(mode="analytical").get_args()
    a3 = get_matuszynska2019(mode="analytical-split").get_args()

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2021() -> None:
    a1 = get_saadat2021(mode="matrix").get_args()
    a2 = get_saadat2021(mode="analytical").get_args()
    a3 = get_saadat2021(mode="analytical-split").get_args()

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


###############################################################################
# High-light tests: matrix and closed-form might diverge at high light intensity
###############################################################################


def test_variants_2016_npq_high_light() -> None:
    a1 = (
        get_matuszynska2016_npq(mode="matrix")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )
    a2 = (
        get_matuszynska2016_npq(mode="analytical")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )
    a3 = (
        get_matuszynska2016_npq(mode="analytical-split")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2016_phd_high_light() -> None:
    a1 = (
        get_matuszynska2016_phd(mode="matrix")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )
    a2 = (
        get_matuszynska2016_phd(mode="analytical")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )
    a3 = (
        get_matuszynska2016_phd(mode="analytical-split")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2019_high_light() -> None:
    a1 = get_matuszynska2019(mode="matrix").update_parameter(n.pfd(), 5000).get_args()
    a2 = (
        get_matuszynska2019(mode="analytical")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )
    a3 = (
        get_matuszynska2019(mode="analytical-split")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)


def test_variants_2021_high_light() -> None:
    a1 = get_saadat2021(mode="matrix").update_parameter(n.pfd(), 5000).get_args()
    a2 = get_saadat2021(mode="analytical").update_parameter(n.pfd(), 5000).get_args()
    a3 = (
        get_saadat2021(mode="analytical-split")
        .update_parameter(n.pfd(), 5000)
        .get_args()
    )

    pd.testing.assert_series_equal(a1, a2[a1.index], atol=1e-12, rtol=1e-12)
    pd.testing.assert_series_equal(a1, a3[a1.index], atol=1e-12, rtol=1e-12)

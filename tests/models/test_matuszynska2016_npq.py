# Matuszyńska 2016 (NPQ)
# https://doi.org/10.1016/j.bbabio.2016.09.003


import pandas as pd
from mxlpy import Simulator

from mxlbricks import get_matuszynska2016npq as get_model


def test_steady_state() -> None:
    model = get_model()
    res = Simulator(model).simulate(100).get_result().unwrap_or_err()
    assert res is not None

    pd.testing.assert_series_equal(
        pd.Series(model.get_initial_conditions()),
        pd.Series(res.get_new_y0()),
    )

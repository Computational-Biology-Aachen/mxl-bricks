from functools import partial

import pandas as pd
from mxlpy import Scipy, Simulator

from mxlbricks import get_yokota1985 as get_model


def test_steady_state() -> None:
    model = get_model()
    res = (
        Simulator(model, integrator=partial(Scipy, method="Radau"))
        .simulate(100)
        .get_result()
        .unwrap_or_err()
    )
    assert res is not None

    pd.testing.assert_series_equal(
        pd.Series(model.get_initial_conditions()),
        pd.Series(res.get_new_y0()),
        rtol=1e-2,
    )

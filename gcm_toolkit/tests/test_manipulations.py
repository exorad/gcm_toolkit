"""testing manipulations functions"""
from gcm_toolkit import GCMT
from gcm_toolkit.tests.test_gcmtools_common import (
    all_raw_testdata,
    all_nc_testdata,
)


def test_horizontal_average(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)
    dsi = tools.get_models()

    area_key = expected.get("area_key", "area_c")
    avg = tools.add_horizontal_average(
        "T", var_key_out="T_g", area_key=area_key
    )

    assert hasattr(dsi, "T_g")
    assert (avg == dsi.T_g).all()
    assert set(dsi.T_g.dims) == {"Z", "time"}


def test_horizontal_overturning(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)
    dsi = tools.get_models()

    v_data = expected.get("v_data", "V")
    psi = tools.add_meridional_overturning(v_data, var_key_out="psi")

    assert hasattr(dsi, "psi")
    assert (psi == dsi.psi).all()
    assert set(dsi.psi.dims) == {"Z", "time", "lat", "lon"}

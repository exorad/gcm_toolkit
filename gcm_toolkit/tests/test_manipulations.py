"""testing manipulations functions"""
from gcm_toolkit import GCMT
from gcm_toolkit.tests.test_gcmtools_common import (
    all_raw_testdata,
    all_nc_testdata,
)


def test_horizontal_average(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on the horizontal average.
    """

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


def test_total_energy(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on the total energy.
    """
    dirname, expected = all_nc_testdata

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname)

    dsi = tools.get_models()

    area_key = expected.get("area_key", "area_c")
    E = tools.add_total_energy(
        var_key_out="E_g", area_key=area_key, temp_key="T"
    )

    assert hasattr(dsi, "E_g")
    assert (E == dsi.E_g).all()
    assert set(dsi.E_g.dims) == {"time"}


def test_total_momentum(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on the total momentum.
    """
    dirname, expected = all_nc_testdata

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname)

    dsi = tools.get_models()

    area_key = expected.get("area_key", "area_c")
    TM = tools.add_total_momentum(
        var_key_out="m_g", area_key=area_key, temp_key="T"
    )

    assert hasattr(dsi, "m_g")
    assert (TM == dsi.m_g).all()
    assert set(dsi.m_g.dims) == {"time"}


def test_rcb(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on the rcb location.
    """
    dirname, expected = all_nc_testdata

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname)

    dsi = tools.get_models()

    area_key = expected.get("area_key", "area_c")
    rcb = tools.add_rcb(
        tol=0.01, var_key_out="rcb", area_key=area_key, temp_key="T"
    )

    assert hasattr(dsi, "rcb")
    assert (rcb == dsi.rcb).all()
    assert set(dsi.rcb.dims) == {"time"}


def test_theta(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on the calculation of theta.
    """
    dirname, expected = all_nc_testdata

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname)

    dsi = tools.get_models()

    theta = tools.add_theta(var_key_out="theta", temp_key="T")

    assert hasattr(dsi, "theta")
    assert (theta == dsi.theta).all()
    assert set(dsi.theta.dims) == {"Z", "time", "lat", "lon"}


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

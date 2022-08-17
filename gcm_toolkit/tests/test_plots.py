"""
General GCMT tests
"""
import os

import pytest
from matplotlib.testing.decorators import image_comparison
from gcm_toolkit import GCMT
from gcm_toolkit.tests.test_gcmtools_common import (
    all_nc_testdata,
    exorad_testdata_nc,
)


@pytest.mark.parametrize("contourf", [True, False])
def test_plot_gcmt_zonal_mean(all_nc_testdata, contourf):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.zonal_mean("U", contourf=contourf)


def test_plot_gcmt_time_evol(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    with pytest.raises(ValueError):
        tools.time_evol("T")

    tools.add_horizontal_average("T", "T_g")
    tools.time_evol("T_g")


@pytest.mark.parametrize("contourf", [True, False])
@pytest.mark.parametrize("lookup_method", ["nearest", "exact", "interpolate"])
def test_plot_gcmt_isobaric_slice(all_nc_testdata, contourf, lookup_method):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.isobaric_slice(
        "T",
        pres=expected["p_domain"][-1],
        lookup_method=lookup_method,
        contourf=contourf,
    )

    with pytest.raises(ValueError):
        tools.isobaric_slice(
            "T",
            pres=expected["p_domain"][-1],
            lookup_method="wrong",
            contourf=contourf,
        )


@image_comparison(["time_evol"], remove_text=True, extensions=[".png"])
def test_plot_gcmt_time_evol_compare(exorad_testdata_nc):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = exorad_testdata_nc

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.add_horizontal_average("T", "T_g")
    tools.time_evol("T_g")


@image_comparison(["isobaric_slice"], remove_text=True, extensions=[".png"])
def test_plot_gcmt_isobaric_slice_compare(
    exorad_testdata_nc, contourf=False, lookup_method="exact"
):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = exorad_testdata_nc

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.isobaric_slice(
        "T",
        pres=expected["p_domain"][-1],
        lookup_method=lookup_method,
        contourf=contourf,
    )


@image_comparison(["zonal_mean"], remove_text=True, extensions=[".png"])
def test_plot_gcmt_zonal_mean_compare(exorad_testdata_nc, contourf=False):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = exorad_testdata_nc

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.zonal_mean("U", contourf=contourf)

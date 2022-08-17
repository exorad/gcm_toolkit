"""
General gcm_toolkit utilities
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.utils.passport import is_the_data_basic, is_the_data_cloudy
from gcm_toolkit.core.units import convert_time, convert_pressure
from gcm_toolkit.tests.test_gcmtools_common import (
    all_raw_testdata,
    all_nc_testdata,
)
import shutil


def test_data_testing_utils(all_nc_testdata):
    dirname, expected = all_nc_testdata
    tag = "readin"

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname, tag=tag)

    ds = tools.models

    assert is_the_data_basic(ds)
    assert not is_the_data_cloudy(ds)

    ds_new = xarray.Dataset()
    ds_new.attrs.update({"tag": tag})

    assert not is_the_data_cloudy(ds_new)
    assert not is_the_data_basic(ds_new)


def test_convert_pressure_failures():
    with pytest.raises(ValueError):
        convert_pressure(xarray.Dataset(), "wrong", "wrong2")


def test_convert_time_failures():
    with pytest.raises(ValueError):
        convert_time(xarray.Dataset(), "wrong", "wrong2")

    with pytest.raises(ValueError):
        convert_time(xarray.Dataset(), "iter", "day")

"""
General gcm_toolkit utilities
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.utils.passport import is_the_data_basic, is_the_data_cloudy
from gcm_toolkit.tests.test_gcmtools_common import all_raw_testdata
import shutil


def test_is_data_testing_utils(all_raw_testdata):
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tag = "raw_readin"

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag=tag)

    ds = tools.models

    assert is_the_data_basic(ds)
    assert not is_the_data_cloudy(ds)

    ds_new = xarray.Dataset()
    ds_new.attrs.update({"tag": tag})

    assert not is_the_data_cloudy(ds_new)
    assert not is_the_data_basic(ds_new)

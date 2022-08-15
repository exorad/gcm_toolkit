"""
Things related to exorad that are not tested in test_gcmt
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import exorad_testdata


def test_exorad(exorad_testdata):
    """Test readin"""

    dirname, expected = exorad_testdata

    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT()
    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, iters="all", tag="0"
    )

    tools.read_raw(
        gcm=expected["gcm"],
        data_path=data_path,
        exclude_iters=expected["iters"],
        tag="empty",
    )
    assert tools.get("empty") is None

    tools.read_raw(
        gcm=expected["gcm"],
        data_path=data_path,
        exclude_iters=expected["iters"][0],
        tag="empty",
    )
    assert tools.get("empty") is None

    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, load_existing=True, tag="0"
    )
    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, load_existing=True, tag="new"
    )

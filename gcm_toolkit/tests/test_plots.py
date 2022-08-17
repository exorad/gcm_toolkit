"""
General GCMT tests
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import all_nc_testdata
from gcm_toolkit.core.units import convert_time
import shutil


def test_plot_gcmt_zonal_mean(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(data_path=dirname)

    tools.get_models()

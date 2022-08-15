"""
Things related to exorad that are not tested in test_gcmt
"""
import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import all_raw_testdata


def test_exorad(all_raw_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

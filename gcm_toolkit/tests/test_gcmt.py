"""
General GCMT tests
"""
import numpy as np
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import all_raw_testdata


def test_create_gcmt_minimal(all_raw_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT()
    tools.read_raw(gcm=expected["gcm"], data_path=data_path)

    dsi = tools.get_models()
    assert isinstance(dsi, xarray.Dataset)


def test_gcmt_get_model_multiple(all_raw_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT()
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag="test1")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag="test2")

    models = tools.get_models()
    assert isinstance(models, GCMDatasetCollection)
    assert models == tools.models

    test1 = tools.get_models("test1")
    test1_one_model = tools.get_one_model("test1")
    assert isinstance(test1, xarray.Dataset)
    assert test1_one_model == test1


def test_gcmt_get_model_single(all_raw_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT()
    tools.read_raw(gcm=expected["gcm"], data_path=data_path)

    models_dict = tools.get_models(always_dict=True)
    assert isinstance(models_dict, GCMDatasetCollection)
    assert set(models_dict.keys()) == {"0"}

    models = tools.get_models()
    assert isinstance(models, xarray.Dataset)
    assert models == tools.models


def test_units(all_raw_testdata):
    """Do tests for the conversion of units in GCMT"""

    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    p_unit = "bar"
    time_unit = "day"

    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_raw(gcm=expected["gcm"], data_path=data_path)

    dsi = tools.get_models()
    assert hasattr(dsi, "time_unit")
    assert hasattr(dsi, "p_unit")
    assert tools.p_unit == p_unit
    assert tools.time_unit == time_unit
    assert dsi.p_unit == p_unit
    assert dsi.time_unit == time_unit

    if p_domain := expected.get("p_domain"):
        # If available in metadata, also check if pressure domain is correct
        p_max = max(p_domain)
        p_min = min(p_domain)
        assert np.isclose(p_max, dsi.Z.max())
        assert np.isclose(p_min, dsi.Z.min())

    if times := expected.get("times"):
        # If available in metadata, also check if timestamps are correct
        assert np.all(np.isclose(dsi.time.values, times))

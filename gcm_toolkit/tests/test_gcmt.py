"""
General GCMT tests
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import all_raw_testdata
import shutil


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
    _ = tools.get_models("test2")
    test1_one_model = tools.get_one_model("test1")
    assert isinstance(test1, xarray.Dataset)
    assert test1_one_model == test1

    assert len(tools) == 2
    assert bool(tools)

    with pytest.raises(ValueError):
        tools.get_one_model()

    with pytest.raises(ValueError):
        _ = tools["wrong"]

    with pytest.raises(ValueError):
        _ = tools[0]

    with pytest.raises(ValueError):
        _ = tools.get(2)

    assert tools.get("test1") == test1
    assert tools.get("wrong", 1) == 1
    assert tools.get("wrong") is None

    get_model_test = list(tools.get_models(always_dict=True).values())
    iter_test = [m for t, m in tools]
    assert get_model_test == iter_test

    _ = tools["test1"]


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

    with pytest.raises(KeyError):
        tools.get_one_model(tag="wrong")

    with pytest.raises(ValueError):
        tools.get_models(tag=2)


def test_units(all_raw_testdata):
    """Do tests for the conversion of units in GCMT"""

    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    p_unit = "bar"
    time_unit = "day"

    with pytest.raises(ValueError):
        _ = GCMT(p_unit="wrong")
    with pytest.raises(ValueError):
        _ = GCMT(time_unit="wrong")

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


def test_set_model(all_raw_testdata):
    """Test the function to set data"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag="raw_readin")

    ds_raw = tools.get("raw_readin")

    with pytest.raises(ValueError):
        tools[2] = ds_raw

    with pytest.raises(ValueError):
        tools["wrong_input"] = 2

    with pytest.raises(ValueError):
        tools["wrong_input"] = ds_raw.drop("T")

    with pytest.raises(NotImplementedError):
        ds_raw.attrs.update({"time_unit": "wrong"})
        tools["wrong_time"] = ds_raw

    ds_raw.attrs.update({"time_unit": tools.time_unit})

    with pytest.raises(NotImplementedError):
        ds_raw.attrs.update({"p_unit": "wrong"})
        tools["wrong_p"] = ds_raw

    ds_raw.attrs.update({"p_unit": tools.p_unit})


def test_load_save(all_raw_testdata):
    """Test load and save"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tag = "raw_readin"

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag=tag)
    ds = tools[tag]
    tools.save(".", method="nc")
    tools.save(".", method="zarr")
    assert os.path.exists(f"{tag}.nc")
    assert os.path.exists(f"{tag}.zarr")

    del tools

    tools = GCMT(write="off")
    tools.load(".", method="nc")
    tools.load(".", method="nc")  # test twice
    os.remove(f"{tag}.nc")

    assert isinstance(tools[tag], xarray.Dataset)
    assert tools[tag] == ds

    del tools

    tools = GCMT(write="off")
    tools.load(".", method="zarr")
    shutil.rmtree(f"{tag}.zarr")

    assert isinstance(tools[tag], xarray.Dataset)
    assert tools[tag] == ds


def test_read_reduced(all_raw_testdata):
    """Test load and save"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tag = "raw_readin"

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag=tag)
    ds = tools[tag]
    tools.save(".", method="nc")
    assert os.path.exists(f"{tag}.nc")

    del tools

    tools = GCMT(write="off")
    tools.read_reduced(f"{tag}.nc", tag=tag)

    assert isinstance(tools[tag], xarray.Dataset)
    assert all(tools[tag].load() == ds)

    os.remove(f"{tag}.nc")


def test_read_incorrect_files(all_raw_testdata):
    """Test load and save"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tag = "raw_readin"

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag=tag)
    ds = tools[tag]
    ds.drop("T").to_netcdf(f"{tag}.nc")
    assert os.path.exists(f"{tag}.nc")

    del tools

    tools = GCMT(write="off")

    with pytest.raises(ValueError):
        tools.read_reduced(f"{tag}.nc", tag=tag)

    os.remove(f"{tag}.nc")

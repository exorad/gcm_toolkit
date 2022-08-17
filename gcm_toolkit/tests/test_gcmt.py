"""
General GCMT tests
"""
import os

import numpy as np
import pytest
import xarray
from gcm_toolkit import GCMT
from gcm_toolkit.gcm_dataset_collection import GCMDatasetCollection
from gcm_toolkit.tests.test_gcmtools_common import (
    all_raw_testdata,
    all_nc_testdata,
)
from gcm_toolkit.core.units import convert_time
import shutil


def test_create_gcmt_minimal(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_nc_testdata
    path = os.path.dirname(dirname)

    tools = GCMT()
    tools.load(path, tag=expected["tag"])

    dsi = tools.get_models()
    assert isinstance(dsi, xarray.Dataset)


def test_gcmt_get_model_multiple(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_nc_testdata

    tools = GCMT()
    tools.read_reduced(dirname, tag="test1")
    tools.read_reduced(dirname, tag="test2")

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


def test_gcmt_get_model_single(all_nc_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""

    dirname, expected = all_nc_testdata
    tools = GCMT()
    tools.read_reduced(dirname)

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


def test_units_reduced(all_nc_testdata):
    """Do tests for the conversion of units in GCMT"""
    dirname, expected = all_nc_testdata
    p_unit = "bar"
    time_unit = "day"

    testfile = "testfile.nc"

    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_reduced(data_path=dirname)

    units_common(expected, tools, p_unit, time_unit)

    dsi = tools.get_models()
    del dsi.attrs["p_unit"]
    dsi.to_netcdf(testfile)

    with pytest.raises(ValueError):
        tools.read_reduced(data_path=testfile)

    del tools
    os.remove(testfile)

    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_reduced(data_path=dirname)
    dsi = tools.get_models()

    del dsi.attrs["time_unit"]
    dsi.to_netcdf(testfile)
    with pytest.raises(ValueError):
        tools.read_reduced(data_path=testfile)

    del tools
    os.remove(testfile)

    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_reduced(data_path=dirname)
    dsi = tools.get_models().copy()
    del dsi.attrs["time_unit"]
    del dsi.attrs["p_unit"]
    dsi.to_netcdf(testfile)

    del tools
    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_reduced(
        data_path=testfile, p_unit_in=p_unit, time_unit_in=time_unit
    )
    units_common(expected, tools, p_unit, time_unit)
    os.remove(testfile)


def test_units_raw(all_raw_testdata):
    """Do tests for the conversion of units in GCMT"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)
    p_unit = "bar"
    time_unit = "day"

    tools = GCMT(p_unit=p_unit, time_unit=time_unit)
    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, load_existing=True
    )

    units_common(expected, tools, p_unit, time_unit)


def units_common(expected, tools, p_unit, time_unit):
    with pytest.raises(ValueError):
        _ = GCMT(p_unit="wrong")
    with pytest.raises(ValueError):
        _ = GCMT(time_unit="wrong")

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
    dsi_new = convert_time(dsi, "day", "iter")
    assert dsi_new.time_unit == "iter"


def test_set_model(all_nc_testdata):
    """Test the function to set data"""
    dirname, expected = all_nc_testdata
    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname, tag="raw_readin")

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


def test_load_save(all_nc_testdata):
    """Test load and save"""
    dirname, expected = all_nc_testdata
    tag = "raw_readin"
    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname, tag=tag)

    ds = tools[tag]
    tools.save(".", method="nc")
    tools.save(".", method="nc", tag=tag)
    tools.save(".", method="zarr")
    tools.save(".", method="zarr", tag=tag)
    tools.save(".", method="zarr", tag=tag, update_along_time=True)

    with pytest.raises(NotImplementedError):
        tools.save(".", method="wrong")

    assert os.path.exists(f"{tag}.nc")
    assert os.path.exists(f"{tag}.zarr")

    del tools

    tools = GCMT(write="off")
    tools.load(".", method="nc")
    tools.load(".", method="nc")  # test twice

    tools.load(".", method="nc", tag=tag)
    tools.load(
        ".", method="nc", tag="wrong"
    )  # will not raise error, and just warn

    os.remove(f"{tag}.nc")

    assert isinstance(tools[tag], xarray.Dataset)
    assert tools[tag] == ds

    del tools

    tools = GCMT(write="off")
    tools.load(".", method="zarr")
    shutil.rmtree(f"{tag}.zarr")

    assert isinstance(tools[tag], xarray.Dataset)
    assert tools[tag] == ds

    with pytest.raises(NotImplementedError):
        tools.load(".", method="wrong")


def test_read_raw(all_raw_testdata):
    """Test load and save"""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tag = "raw_readin"

    tools = GCMT(write="off")
    tools.read_raw(gcm=expected["gcm"], data_path=data_path, tag=tag)

    assert isinstance(tools[tag], xarray.Dataset)


def test_read_incorrect_files(all_nc_testdata):
    """Test load and save"""
    dirname, expected = all_nc_testdata
    tag = "nc_incorrect"

    tools = GCMT(write="off")
    tools.read_reduced(data_path=dirname, tag=tag)
    ds = tools[tag]
    ds.drop("T").to_netcdf(f"{tag}.nc")
    assert os.path.exists(f"{tag}.nc")

    del tools

    tools = GCMT(write="off")

    with pytest.raises(ValueError):
        tools.read_reduced(f"{tag}.nc", tag=tag)

    os.remove(f"{tag}.nc")


def test_gcmt_read_raw_existing(all_raw_testdata):
    """Create a minimal gcm_toolkit object and do simple tests on it."""
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    tools = GCMT()
    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, load_existing=True
    )
    assert bool(tools), "did not load a dataset!"


def test_read_incorrect_raw_gcm():
    """Raise error, when we try to read a GCM that is not yet supported"""
    tools = GCMT(write="off")

    with pytest.raises(NotImplementedError):
        tools.read_raw("wrong_gcm", data_path=".")

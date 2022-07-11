import numpy as np
import xarray

from GCMtools.GCMDatasetCollection import GCMDatasetCollection
from GCMtools.tests.test_gcmtools_common import *
from GCMtools import GCMT


def test_create_gcmt_minimal(all_raw_testdata):
    """Create a minimal gcmt object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    gcmt = GCMT()
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path)

    ds = gcmt.get_models()
    assert isinstance(ds, xarray.Dataset)


def test_gcmt_get_model_multiple(all_raw_testdata):
    """Create a minimal gcmt object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    gcmt = GCMT()
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path, tag='test1')
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path, tag='test2')

    models = gcmt.get_models()
    assert isinstance(models, GCMDatasetCollection)
    assert models == gcmt.models

    test1 = gcmt.get_models('test1')
    test1_one_model = gcmt.get_one_model('test1')
    assert isinstance(test1, xarray.Dataset)
    assert test1_one_model == test1


def test_gcmt_get_model_single(all_raw_testdata):
    """Create a minimal gcmt object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    gcmt = GCMT()
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path)

    models_dict = gcmt.get_models(always_dict=True)
    assert isinstance(models_dict, GCMDatasetCollection)
    assert set(models_dict.keys()) == {'0'}

    models = gcmt.get_models()
    assert isinstance(models, xarray.Dataset)
    assert models == gcmt.models


def test_units(all_raw_testdata):
    """Do tests for the conversion of units in GCMT"""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    p_unit = 'bar'
    time_unit = 'day'

    gcmt = GCMT(p_unit=p_unit, time_unit=time_unit)
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path)

    ds = gcmt.get_models()
    assert hasattr(ds, 'time_unit')
    assert hasattr(ds, 'p_unit')
    assert gcmt.p_unit == p_unit
    assert gcmt.time_unit == time_unit
    assert ds.p_unit == p_unit
    assert ds.time_unit == time_unit

    if p_domain := expected.get('p_domain'):
        # If available in metadata, also check if pressure domain is correct
        p_max = max(p_domain)
        p_min = min(p_domain)
        assert np.isclose(p_max, ds.Z.max())
        assert np.isclose(p_min, ds.Z.min())

    if times := expected.get('times'):
        # If available in metadata, also check if timestamps are correct
        assert np.all(np.isclose(ds.time.values, times))




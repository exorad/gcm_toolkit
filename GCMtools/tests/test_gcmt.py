import xarray

from GCMtools.tests.test_gcmtools_common import *
from GCMtools import GCMT


def test_create_gcmt_minimal(all_raw_testdata):
    """Create a minimal xarray object with only dimensions in it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    gcmt = GCMT(p_unit='bar', time_unit='day', write='off')
    gcmt.read_raw(gcm=expected['gcm'], data_path=data_path)

    ds = gcmt.get_models()
    assert isinstance(ds, xarray.Dataset)
    assert hasattr(ds, 'time_unit')
    assert hasattr(ds, 'p_unit')




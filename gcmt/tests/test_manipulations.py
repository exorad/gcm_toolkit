from gcmt import GCMT


def test_horizontal_average(all_raw_testdata):
    """Create a minimal gcmt object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    tools = GCMT()
    tools.read_raw(gcm=expected['gcm'], data_path=data_path)
    ds = tools.get_models()

    area_key = expected.get('area_key', 'area_c')
    avg = tools.add_horizontal_average('T', var_key_out='T_g', area_key=area_key)

    assert hasattr(ds, 'T_g')
    assert (avg == ds.T_g).all()
    assert set(ds.T_g.dims) == {'Z', 'time'}


def test_horizontal_overturning(all_raw_testdata):
    """Create a minimal gcmt object and do simple tests on it."""

    dirname, expected = all_raw_testdata
    data_path = expected.get('rel_data_dir', '{}').format(dirname)

    tools = GCMT()
    tools.read_raw(gcm=expected['gcm'], data_path=data_path)
    ds = tools.get_models()

    v_data = expected.get('v_data', 'V')
    psi = tools.add_meridional_overturning(v_data, var_key_out='psi')

    assert hasattr(ds, 'psi')
    assert (psi == ds.psi).all()
    assert set(ds.psi.dims) == {'Z', 'time', 'lat', 'lon'}

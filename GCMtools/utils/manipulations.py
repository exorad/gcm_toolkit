from GCMtools.core.const import VARNAMES as c

def m_add_horizontal_average(gcmt, var_key, var_key_out=None, area_key='area_c', tag=None):
    """
    Calculate horizontal averaged quantities. Horizontal averages
    are calculated as:

        \bar q = \int{q dA}/\int{dA}

    Parameters
    ----------
    var_key: str
        The key of the variable quantity that should be plotted.
    tag : str, optional
        The tag of the dataset that should be used. If no tag is provided,
        and multiple datasets are available, an error is raised.
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the averages and not change the dataset inplace.
    area_key: str, optional
        Variable key in the dataset for the area of grid cells

    Returns
    -------
    TODO
    """
    # print information
    wrt.write_status('STAT', 'Calculate horizontal average')
    wrt.write_status('INFO', 'Tag: ' + tag)
    wrt.write_status('INFO', 'Varaible to be plotted: ' + var_key)
    wrt.write_status('INFO', 'Output variable: ' + var_key_out)
    wrt.write_status('INFO', 'Area of grid cells: ' + area_key)

    ds = gcmt._models.get_one_model(tag)
    avg = (ds[area_key]*ds[var_key]).sum(dim=[c['lon'],c['lat']])/ds[area_key].sum(dim=[c['lon'],c['lat']])

    if var_key_out is not None:
        ds.update({var_key_out: avg})

    return avg
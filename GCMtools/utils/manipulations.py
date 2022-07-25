import numpy as np

import GCMtools.core.writer as wrt
from GCMtools.core.const import VARNAMES as c


def m_add_horizontal_average(ds, var_key, var_key_out=None, area_key='area_c'):
    """
    Calculate horizontal averaged quantities. Horizontal averages
    are calculated as area-weighted averages of quantity q:

        \bar q = \int{q dA}/\int{dA}

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset for which the calculation should be performed
    var_key: str
        The key of the variable quantity that should be plotted.
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the averages and not change the dataset inplace.
    area_key: str, optional
        Variable key in the dataset for the area of grid cells

    Returns
    -------
    avg : xarray.DataArray
        A dataArray with reduced dimensionality, containing the horizontally 
        averaged quantity.
    """
    # print information
    wrt.write_status('STAT', 'Calculate horizontal average')
    wrt.write_status('INFO', 'Variable to be plotted: ' + var_key)
    if var_key_out is not None:
        wrt.write_status('INFO', 'Output variable: ' + var_key_out)
    wrt.write_status('INFO', 'Area of grid cells: ' + area_key)

    avg = (ds[area_key] * ds[var_key]).sum(dim=[c['lon'], c['lat']]) / ds[area_key].sum(dim=[c['lon'], c['lat']])

    if var_key_out is not None:
        ds.update({var_key_out: avg})

    return avg


def m_add_meridional_overturning(ds, v_data='V', var_key_out=None):
    """
    Calculate meridional overturning streamfunction. This quantity psi is 
    computed by integrating the zonal-mean meridional velocity \bar V along
    pressure, and weighting with 2*pi*R_p / g times the cosine of latitude,
    where R_p is the planetary radius and g is the surface gravity:
    
    \bar \psi = 2 \pi R_p / g \cos{lat} \int{\bar V dp'}

    (see e.g. Carone et al. (2018), Eq. 7)

    Parameters
    ----------
    ds: xarray.Dataset
        The dataset for which the calculation should be performed
    v_data: str
        The key that holds the data to meridional velocity
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the overturning circulation and not change the dataset inplace.
    """

    # print information
    wrt.write_status('STAT', 'Calculate meridional overturning streamfunction')
    if var_key_out is not None:
        wrt.write_status('INFO', 'Output variable: ' + var_key_out)

    V_integral = ds[v_data].cumulative_integrate(coord='Z')

    if ds.attrs.get('p_unit') == 'bar':
        # convert to SI, if needed
        V_integral = V_integral / 1.0e5

    psi = 2 * np.pi * np.cos(ds.lat / 180 * np.pi) * ds.R_p / ds.g * V_integral

    if var_key_out is not None:
        ds.update({var_key_out: psi})

    return psi

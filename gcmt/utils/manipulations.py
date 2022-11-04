"""
Functions to manipulate GCM data
"""
import numpy as np
from ..core import writer as wrt
from ..core.const import VARNAMES as c


def m_add_horizontal_average(dsi, var_key, var_key_out=None, area_key='area_c'):
    """
    Calculate horizontal averaged quantities. Horizontal averages
    are calculated as area-weighted averages of quantity q:

        bar q = int{q dA}/int{dA}

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

    avg = ((dsi[area_key] * dsi[var_key]).sum(dim=[c['lon'], c['lat']]) /
           dsi[area_key].sum(dim=[c['lon'], c['lat']]))

    if var_key_out is not None:
        dsi.update({var_key_out: avg})

    return avg


def m_add_meridional_overturning(dsi, v_data='V', var_key_out=None):
    """
    Calculate meridional overturning streamfunction. This quantity psi is
    computed by integrating the zonal-mean meridional velocity \bar V along
    pressure, and weighting with 2*pi*R_p / g times the cosine of latitude,
    where R_p is the planetary radius and g is the surface gravity:

    bar{psi} = 2 pi R_p / g cos{lat} int{bar{V} dp'}

    (see e.g. Carone et al. (2018), Eq. 7)

    Parameters
    ----------
    dsi: xarray.Dataset
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

    v_integral = dsi[v_data].cumulative_integrate(coord='Z')

    if dsi.attrs.get('p_unit') == 'bar':
        # convert to SI, if needed
        v_integral = v_integral / 1.0e5

    psi = 2 * np.pi * np.cos(dsi.lat / 180 * np.pi) * dsi.R_p / dsi.g * v_integral

    if var_key_out is not None:
        dsi.update({var_key_out: psi})

    return psi


def m_extend_upward(temp_array, p_low, method=None, n_p=20, T_therm=None,
                    p_therm_high=None):
    """
    Extend the given temperature DataArray upward to the given pressure value.
    Multiple methods can be used to extend the data upwards:
        - 'isothermal' copies the temperature at lowest pressure and adds layers
          of the same temperature on top
        - 'thermosphere' imposes a simple parametrized thermosphere with a hot
          day-side and gradual decrease to an isothermal night-side
          (cfr. Baeyens+2022, Section 4.2)

    Parameters
    ----------
    temp_array: DataArray
        The temperature data-array that is to be extended to lower pressures.
    p_low: float or list
        If a float is given, extend to the given low-pressure limit using an
        arbitrarily chosen pressure-spacing. If a list is given, use these as
        pressures to compute temperature extension.
    method: string
        Method by which the temperature is to be extended.
    n_p: int, optional
        Number of pressure layers that should be added (only used in case p_low
        is a single number).
    T_therm: float, optional
        Maximal day-side temperature for the upper thermosphere (required for
        the 'thermosphere' method.)
    p_therm_high: float, optional
        The high-pressure limit at which thermospheric heating should start, in
        native pressure units (required for the 'thermosphere' method).

    Returns
    -------
    temp_array_ext: DataArray
        The extended temperature array.
    """
    import xarray as xr

    # Check if the input pressure is suitable
    p_low_current = min(temp_array.Z.values) # lowest pressure in current data
    if hasattr(p_low, "__len__"):           # p_low is a list
        p_low = [ip for ip in p_low if ip < p_low_current]
        p_low = sorted(p_low)[::-1]
    elif isinstance(p_low, (int, float)):   # p_low is a scalar
        p_low = np.logspace(np.log10(p_low_current), np.log10(p_low), n_p+1)
        p_low = p_low[1:]
    else:                                   # p_low is something else
        wrt.write_status('ERROR', 'p_low should be a float, int, list, or array.')

    # Now p_low consists of a sorted list of low pressures

    # Construct an extended dataArray (initially filled with zeros)
    new_coords = dict(time= temp_array.time, Z=p_low, lat=temp_array.lat,
                      lon=temp_array.lon)
    # NB: in order to avoid problems with merging the extended and the original
    # temperature array, we need to copy over all coordinates. We just assign
    # them NaN values since we aren't interested in their values.
    for c in temp_array.coords:
        if c not in new_coords:
            new_coords[c] = np.nan
    nz = len(new_coords['Z'])
    nlat = len(new_coords['lat'])
    nlon = len(new_coords['lon'])
    try:
        nt = len(new_coords['time'])
        time_coord = True # flag for an explicit time coordinate --> 4D array
        extended_values = np.zeros((nt, nz, nlat, nlon))
    except TypeError:
        nt = 1
        time_coord = False # no explicit time coordinate --> 3D array
        extended_values = np.zeros((nz, nlat, nlon))

    # Fill the extended data according to the chosen method
    # ISOTHERMAL
    if method == 'isothermal':
        iso_slice = temp_array.isel(Z=-1).values    # constant upper temperature
        for i in range(0, nz):
            if time_coord:
                extended_values[:,i,:,:] = iso_slice # copy for each new pressure
            else:
                extended_values[i,:,:] = iso_slice # copy for each new pressure

    # THERMOSPHERE
    elif method == 'thermosphere':
        import math
        # Check if all parameters are given
        if T_therm is None:
            wrt.write_status('ERROR', 'Argument "T_therm" is required for the thermosphere method.')
        if p_therm_high is None:
            wrt.write_status('ERROR', 'Argument "p_therm_high" is required for the thermosphere method.')
        p_therm_low = min(p_low)   # minimum pressure for thermosphere extention

        # For each time step...
        for it in range(0, nt):
            # and each vertical column...
            for ilon in range(0, nlon):
                for ilat in range(0, nlat):
                    # ... save the top temperature as minimum
                    T_min = temp_array.isel(Z=-1, lat=ilat, lon=ilon)
                    # ... calculate angle of incidence
                    mu = math.cos(math.radians(T_min.lon)) * \
                         abs( math.cos(math.radians(T_min.lat)) )
                    # ... determine the 'thermospheric' upper temperature.
                    T_max = max(T_min, T_therm * mu)

                    # For each pressure level...
                    for i in range(0, nz):
                        ip = p_low[i]
                        # ... if the pressure is below the thermosphere...
                        if ip < p_therm_high and mu >= 0:
                            # ... do the interpolation to T_max
                            T_new = T_min + (T_max-T_min)/(math.log10(p_therm_high)-math.log10(p_therm_low)) * (math.log10(p_therm_high)-math.log10(ip))
                        else: # ... otherwise don't change anything.
                            T_new = T_min
                        if time_coord:
                            extended_values[it, i, ilat, ilon] = T_new
                        else:
                            extended_values[i, ilat, ilon] = T_new
    else:
        msg = 'Specify either "isothermal" or "thermosphere" as method.'
        wrt.write_status('ERROR', msg)

    # Make the new extended data into a DataArray
    extension = xr.DataArray(extended_values,
                             coords=new_coords,
                             dims=temp_array.dims)

    # Combine the extended and original data
    temp_array = xr.concat([temp_array, extension], dim='Z',
                           compat='no_conflicts')

    return temp_array

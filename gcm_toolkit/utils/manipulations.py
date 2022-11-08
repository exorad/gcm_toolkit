"""
Functions to manipulate GCM data
"""
import numpy as np
import xarray as xr

from ..core import writer as wrt
from ..core.const import VARNAMES as c
from ..core.units import convert_pressure


def m_add_horizontal_average(
    dsi, var_key, var_key_out=None, part="global", area_key="area_c"
):
    """
    Calculate horizontal averaged quantities. Horizontal averages
    are calculated as area-weighted averages of quantity q:

        bar q = int{q dA}/int{dA}

    Parameters
    ----------
    dsi: xarray.Dataset
        The dataset for which the calculation should be performed
    var_key: str, xarray.DataArray
        The key or array of the variable quantity that should be averaged.
        If str, it will try to look up the key in the dataset.
        If DataArray, it will use this one instead.
    var_key_out: str, optional
        variable name used to store the outcome.
        If not provided, this script will just
        return the averages and not change the dataset inplace.
    part: dict or str, optional
        'global': global average
        'night': only nightside (defined around +-180,0)
        'day': only dayside (defined around 0,0)
        'morning': morning terminator (average around lon=[-100,-80])
        'evening': evening terminator (average around lon=[80,100])
        Alternatively you may specify a dict in the following way
        (example for morn. term.):
        part = {'lon': [-100,-80], 'lat':[-90,90], 'inv':False}
        The 'lon', 'lat' specify regions of lon and lat that should be used,
        whereas 'inv' (optional, default False) gives the option to invert
        the lon and lat regions (e.g., exclude instead of include for average)
    area_key: str, optional
        Variable key in the dataset for the area of grid cells

    Returns
    -------
    avg : xarray.DataArray
        A dataArray with reduced dimensionality, containing the horizontally
        averaged quantity.
    """
    # print information
    wrt.write_status("STAT", "Calculate horizontal average")
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)
    wrt.write_status("INFO", "Area of grid cells: " + area_key)

    if isinstance(var_key, str):
        data = dsi[var_key]
        wrt.write_status("INFO", "Variable to be averaged: " + var_key)
    elif isinstance(var_key, xr.DataArray):
        data = var_key
        wrt.write_status("INFO", "Variable to be averaged is taken from input")
    else:
        raise ValueError(
            "var_key needs to be either str (key in Dataset) or DataArray"
        )

    # Determine the area over which we want to average:
    if isinstance(part, str):
        if part == "global":
            part_internal = {}
            wrt.write_status("INFO", "Performing global average")
        elif part == "day":
            part_internal = {"lon": [-90, 90]}
            wrt.write_status("INFO", "Performing dayside average")
        elif part == "night":
            part_internal = {"lon": [-90, 90], "inv": True}
            wrt.write_status("INFO", "Performing nightside average")
        elif part == "evening":
            part_internal = {"lon": [80, 100]}
            wrt.write_status("INFO", "Performing evening terminator average")
        elif part == "morning":
            part_internal = {"lon": [-100, -80]}
            wrt.write_status("INFO", "Performing morning terminator average")
        else:
            raise ValueError(
                "If you specify a string for part, it needs to be either"
                " morning, evening, day, night or global"
            )
    elif isinstance(part, dict):
        wrt.write_status("INFO", "Performing user specified average")
        part_internal = part
    else:
        raise ValueError("Please use a dict or a string for part.")

    if lon := part_internal.get("lon"):
        lon_bool = np.logical_and(
            dsi[c["lon"]] <= max(lon),
            dsi[c["lon"]] >= min(lon),
        )
    else:
        lon_bool = xr.ones_like(dsi[area_key], dtype=np.bool)

    if lat := part_internal.get("lat"):
        lat_bool = np.logical_and(
            dsi[c["lat"]] <= max(lat),
            dsi[c["lat"]] >= min(lat),
        )
    else:
        lat_bool = xr.ones_like(dsi[area_key], dtype=np.bool)

    combined_bool = np.logical_and(lat_bool, lon_bool)

    if part_internal.get("inv", False):
        area = xr.where(combined_bool, 0, dsi[area_key])
    else:
        area = xr.where(combined_bool, dsi[area_key], 0)

    avg = (area * data).sum(dim=[c["lon"], c["lat"]]) / area.sum(
        dim=[c["lon"], c["lat"]]
    )

    if var_key_out is not None:
        dsi.update({var_key_out: avg})

    return avg


def m_add_total_energy(
    dsi, var_key_out=None, area_key="area_c", temp_key="T", return_all=False
):
    """
    Calculate the total Energy of the GCM. See e.g.,
    https://ui.adsabs.harvard.edu/abs/2014Icar..229..355P, Eq. 16

    Parameters
    ----------
    dsi: xarray.Dataset
        The dataset for which the calculation should be performed
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the averages and not change the dataset inplace.
    area_key: str, optional
        Variable key in the dataset for the area of grid cells
    temp_key: str, optional
        The key to look up the temperature
    return_all: bool, optional
        Also return the partial energies

    Returns
    -------
    tot_energy : xarray.DataArray
        A dataArray with reduced dimensionality, containing the total energy.
    therm_energy : xarray.DataArray, optional
        A dataArray with reduced dimensionality, containing the thermal energy.
    pot_energy : xarray.DataArray, optional
        A dataArray with reduced dimensionality, containing the potential energy.
    kin_energy : xarray.DataArray, optional
        A dataArray with reduced dimensionality, containing the kinetic energy.

    """
    # print information
    wrt.write_status("STAT", "Calculate total energy")
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)
    wrt.write_status("INFO", "Area of grid cells: " + area_key)
    wrt.write_status("INFO", "Temperature variable: " + temp_key)

    dsi_calc = dsi.copy()
    dsi_calc = convert_pressure(dsi_calc, dsi_calc.p_unit, "Pa")

    dzdp, rho = _calc_hydrostat_eq(dsi_calc, temp_key)
    z_geo = dzdp.cumulative_integrate(coord=c["Z"])

    d_pot_energy = z_geo * dsi_calc.attrs[c["g"]]

    d_kin_energy = 0.5 * (dsi_calc[c["U"]] ** 2 + dsi_calc[c["V"]] ** 2)

    d_therm_energy = dsi_calc.attrs[c["cp"]] * dsi_calc[temp_key]

    kin_energy = _integrate_over_mass(
        quant_to_int=d_kin_energy, area=dsi_calc[area_key], dzdp=dzdp, rho=rho
    )
    therm_energy = _integrate_over_mass(
        quant_to_int=d_therm_energy,
        area=dsi_calc[area_key],
        dzdp=dzdp,
        rho=rho,
    )
    pot_energy = _integrate_over_mass(
        quant_to_int=d_pot_energy, area=dsi_calc[area_key], dzdp=dzdp, rho=rho
    )

    tot_energy = kin_energy + therm_energy + pot_energy

    if var_key_out is not None:
        dsi.update({var_key_out: tot_energy})
        dsi.update({var_key_out + "_th": therm_energy})
        dsi.update({var_key_out + "_pot": pot_energy})
        dsi.update({var_key_out + "_kin": kin_energy})

    if return_all:
        return tot_energy, therm_energy, pot_energy, kin_energy
    else:
        return tot_energy


def _integrate_over_mass(quant_to_int, area, dzdp, rho):
    """Helper function that carries out a mass integral (dM = rho dV)."""
    return (
        (quant_to_int * rho * area * dzdp)
        .sum(dim=[c["lon"], c["lat"]])
        .integrate(coord=c["Z"])
    )


def _calc_hydrostat_eq(dsi, temp_key):
    """
    Helper function that calculated rho and dz/dp
    based on the ideal gas equation and on hydrostatic eq.
    """
    rho = dsi[c["Z"]] / dsi.attrs[c["R"]] / dsi[temp_key]
    dzdp = -1 / rho / dsi.attrs[c["g"]]
    return dzdp, rho


def m_add_total_momentum(
    dsi, var_key_out=None, area_key="area_c", temp_key="T"
):
    """
    Calculate the total angular momentum of the GCM. See e.g.,
    https://ui.adsabs.harvard.edu/abs/2014Icar..229..355P, Eq. 17

    Parameters
    ----------
    dsi: xarray.Dataset
        The dataset for which the calculation should be performed
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the averages and not change the dataset inplace.
    area_key: str, optional
        Variable key in the dataset for the area of grid cells
    temp_key: str, optional
        The key to look up the temperature (needed for density calculation)

    Returns
    -------
    momentum : xarray.DataArray
        A dataArray with reduced dimensionality, containing the total momentum.
    """
    # print information
    wrt.write_status("STAT", "Calculate total angular momentum")
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)
    wrt.write_status("INFO", "Area of grid cells: " + area_key)
    wrt.write_status("INFO", "Temperature variable: " + temp_key)

    dsi_calc = dsi.copy()
    dsi_calc = convert_pressure(dsi_calc, dsi_calc.p_unit, "Pa")

    dzdp, rho = _calc_hydrostat_eq(dsi_calc, temp_key)
    cosphi = np.cos(dsi_calc[c["lat"]] * np.pi / 180)
    d_momentum = (
        (
            2
            * np.pi
            / dsi_calc.attrs[c["P_rot"]]
            * dsi_calc.attrs[c["R_p"]]
            * cosphi
            + dsi_calc[c["U"]]
        )
        * cosphi
        * dsi_calc.attrs[c["R_p"]]
    )
    momentum = _integrate_over_mass(
        quant_to_int=d_momentum, area=dsi_calc[area_key], dzdp=dzdp, rho=rho
    )

    if var_key_out is not None:
        dsi.update({var_key_out: momentum})

    return momentum


def m_add_theta(dsi, var_key_out=None, temp_key="T"):
    """
    Convert temperature to potential temperature with respect to model boundary.

    Parameters
    ----------
    dsi: xarray.Dataset
        The dataset for which the calculation should be performed
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return theta and not change the dataset inplace.
    temp_key: str, optional
        The key to look up the temperature

    Returns
    -------
    theta : xarray.DataArray
        A dataArray with reduced dimensionality, containing the potential temperature
    """
    theta = dsi[c[temp_key]] * (dsi[c["Z"]].max() / dsi[c["Z"]]) ** (
        dsi.attrs[c["R"]] / dsi.attrs[c["cp"]]
    )

    if var_key_out is not None:
        dsi.update({var_key_out: theta})

    return theta


def m_add_rcb(
    dsi,
    tol=0.01,
    var_key_out=None,
    part="global",
    area_key="area_c",
    temp_key="T",
):
    """
    Calculate the radiative convective boundary (rcb) by searching
    (from the bottom upwards) for the first occurance of a deviation
    from an adiabatic temperature profile.

    Operates on and calculates the horizontal average.

    Parameters
    ----------
    dsi: xarray.Dataset
        The dataset for which the calculation should be performed
    tol: float
        tolerance for the relative deviation from adiabat
    var_key_out: str, optional
        variable name used to store the outcome. If not provided, this script will just
        return the averages and not change the dataset inplace.
    part: dict or str, optional
        'global': global average
        'night': only nightside (defined around +-180,0)
        'day': only dayside (defined around 0,0)
        'morning': morning terminator (average around lon=[-100,-80])
        'evening': evening terminator (average around lon=[80,100])
        Alternatively you may specify a dict in the following way
        (example for morn. term.):
        part = {'lon': [-100,-80], 'lat':[-90,90], 'inv':False}
        The 'lon', 'lat' specify regions of lon and lat that should be used,
        whereas 'inv' (optional, default False) gives the option to invert
        the lon and lat regions (e.g., exclude instead of include for average)
    area_key: str, optional
        Variable key in the dataset for the area of grid cells
    temp_key: str, optional
        The key to look up the temperature (needed for density calculation)

    Returns
    -------
    rcb : xarray.DataArray
        A dataArray with reduced dimensionality,
        containing the pressure of the rcb location.
    """
    # print information
    wrt.write_status("STAT", "Calculate the location of the rcb")
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)
    wrt.write_status("INFO", "Area of grid cells: " + area_key)
    wrt.write_status("INFO", "Temperature variable: " + temp_key)

    dsi_calc = dsi.copy()
    dsi_calc = convert_pressure(dsi_calc, dsi_calc.p_unit, "Pa")

    m_add_theta(dsi_calc, temp_key=temp_key, var_key_out="theta")
    theta_g = m_add_horizontal_average(dsi_calc, var_key="theta", part=part)

    rcb_loc = (
        abs(
            (theta_g - theta_g.isel(**{c["Z"]: 0}))
            / theta_g.isel(**{c["Z"]: 0})
        )
        < tol
    ).argmin(dim=c["Z"])

    rcb = dsi[c["Z"]].isel(**{c["Z"]: rcb_loc})

    if var_key_out is not None:
        dsi.update({var_key_out: rcb})

    return rcb


def m_add_meridional_overturning(dsi, v_data="V", var_key_out=None):
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
    wrt.write_status("STAT", "Calculate meridional overturning streamfunction")
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)

    v_integral = dsi[v_data].cumulative_integrate(coord="Z")

    if dsi.attrs.get("p_unit") == "bar":
        # convert to SI, if needed
        v_integral = v_integral / 1.0e5

    psi = (
        2
        * np.pi
        * np.cos(dsi.lat / 180 * np.pi)
        * dsi.R_p
        / dsi.g
        * v_integral
    )

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

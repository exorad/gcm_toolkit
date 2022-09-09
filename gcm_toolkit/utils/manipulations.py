"""
Functions to manipulate GCM data
"""
import numpy as np
from ..core import writer as wrt
from ..core.const import VARNAMES as c
from ..core.units import convert_pressure


def m_add_horizontal_average(
    dsi, var_key, var_key_out=None, area_key="area_c"
):
    """
    Calculate horizontal averaged quantities. Horizontal averages
    are calculated as area-weighted averages of quantity q:

        bar q = int{q dA}/int{dA}

    Parameters
    ----------
    dsi: xarray.Dataset
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
    wrt.write_status("STAT", "Calculate horizontal average")
    wrt.write_status("INFO", "Variable to be plotted: " + var_key)
    if var_key_out is not None:
        wrt.write_status("INFO", "Output variable: " + var_key_out)
    wrt.write_status("INFO", "Area of grid cells: " + area_key)

    avg = (dsi[area_key] * dsi[var_key]).sum(dim=[c["lon"], c["lat"]]) / dsi[
        area_key
    ].sum(dim=[c["lon"], c["lat"]])

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
    dsi, tol=0.01, var_key_out=None, area_key="area_c", temp_key="T"
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
    theta_g = m_add_horizontal_average(dsi_calc, var_key="theta")

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

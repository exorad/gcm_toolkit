"""
utils to work with cubedsphere
"""
import os

import cubedsphere as cs
import cubedsphere.const as c
import numpy as np
import xgcm
from f90nml import Parser
from xgcm.autogenerate import generate_grid_ds


class MITgcmDataParser(Parser):
    """
    MIT gcm data parser class to remember common variables
    """

    def __init__(self):
        super().__init__()
        self.comment_tokens += "#"
        self.end_comma = True
        self.indent = " "
        self.column_width = 72
        self.sparse_arrays = True


def get_parameter(datafile, keyword, default=None):
    """
    Function to parse the MITgcm 'data' file and return the parameter values
    of the given specific keyword.

    Parameters
    ----------
    datafile: string
        Full path to the MITgcm data file.
    keyword: string
        Parameter of which the value is required.

    Returns
    ----------
    value: string
        The value associated with the given keyword is returned as a string (!).
    """

    if not os.path.isfile(datafile):
        raise FileNotFoundError("could not find the datafile.")

    parser = MITgcmDataParser()
    data = parser.read(datafile)

    for section in data:
        for key, val in data[section].items():
            if key.lower() == keyword.lower():
                return val

    if default is None:
        raise KeyError(f"Keyword {keyword} not found in datafile.")
    return default


def convert_winds_and_t(dsi, temp_dim, w_dim):
    """
    Convert winds and temperature in dataset.
    Winds are converted from Pa/s to m/s.
    Temperatures are converted from potential temperature to ordinary temperature

    Parameters
    ----------
    dsi: Dataset
        dataset to be converted
    temp_dim: str
        temperature datadimension to be converted
    w_dim: str
        vertical wind datadimension to be converted

    Returns
    -------
    ds: Dataset
        dataset with converted dimension
    """
    kappa = dsi.attrs["R"] / dsi.attrs["cp"]
    dsi[temp_dim] = dsi[temp_dim] * (dsi[c.Z] / dsi.attrs["p_ref"]) ** kappa

    # calculate scale height
    h_val = dsi.attrs["R"] / dsi.attrs["g"] * dsi[temp_dim]

    # calculate geometric height
    dsi[c.Z_geo] = -h_val * np.log(dsi[c.Z] / dsi.attrs["p_ref"])

    if w_dim in dsi:
        # interpolate vertical windspeed to cell center:
        if c.FACEDIM in dsi.dims:
            grid = cs.init_grid_CS(ds=dsi)
        else:
            grid = cs.init_grid_LL(ds=dsi)

        w_interp = grid.interp(dsi[w_dim], axis=c.Z, to="center")

        # convert vertical wind speed from Pa/s to m/s
        dsi[w_dim] = -w_interp * h_val / dsi[c.Z]

    return dsi


def exorad_postprocessing(dsi, outdir=None, datafile=None):
    """
    Preliminaray postprocessing on exorad dataset.
    This function converts the vertical windspeed from Pa
    into meters and saves attributes to the dataset.

    Parameters
    ----------
    dsi: Dataset
        dataset to be extended
    outdir: string
        directory in which to find the data file
        (following the convention f'{outdir}/data')
    datafile: string
        alternatively specify datafile directly

    Returns
    ----------
    ds:
        Dataset to be returned
    """
    stri = "please specify a datafile or a folder where we can find a datafile"
    assert outdir is not None or datafile is not None, stri

    if outdir is not None:
        datafile = f"{outdir}/data"

    # Add metadata
    radius = float(
        get_parameter(datafile, "rSphere", 6370e3)
    )  # planet radius in m
    p_rot = float(
        get_parameter(datafile, "rotationperiod", 8.6164e4)
    )  # planet rotationperiod in m
    attrs = {
        "p_ref": float(
            get_parameter(datafile, "Ro_SeaLevel", 1.0e5)
        ),  # bot layer pres in pa
        "cp": float(
            get_parameter(datafile, "atm_Cp", 1.004e3)
        ),  # heat cap at constant pres
        "R": float(
            get_parameter(datafile, "atm_Rd", 2.868571e2)
        ),  # specific gas constant
        "g": float(
            get_parameter(datafile, "gravity", 9.81)
        ),  # surface gravity in m/s^2
        "dt": int(
            get_parameter(datafile, "deltaT", 0.0)
        ),  # time step size in s
        "R_p": radius,
        "P_rot": p_rot,
        "P_orb": p_rot,  # change, when this is available
    }

    dsi.attrs.update(attrs)

    # Convert Temperature and winds
    if c.T in dsi:
        dsi = convert_winds_and_t(dsi, c.T, c.W)
    if c.Ttave in dsi:
        dsi = convert_winds_and_t(dsi, c.Ttave, c.wVeltave)

    # Add metrics to dataset
    if c.FACEDIM not in dsi.dims:
        dsi = add_distances(dsi, radius=radius)

    return dsi


def add_distances(dsi, radius):
    """Add metric distances into dataset if dataset is in lon lat.

    PARAMETERS
    ----------
    ds : xarray.DataSet (in lon,lat)
    radius: planetary radius im meters

    RETURNS
    -------
    ds  : xarray.DataArray distance inferred from dlon
    dy  : xarray.DataArray distance inferred from dlat
    """

    def dll_dist(dlon, dlat, lon, lat, radius):
        """Converts lat/lon differentials into distances in meters

        PARAMETERS
        ----------
        dlon : xarray.DataArray longitude differentials
        dlat : xarray.DataArray latitude differentials
        lon  : xarray.DataArray longitude values
        lat  : xarray.DataArray latitude values
        radius: planetary radius im meters

        RETURNS
        -------
        dx  : xarray.DataArray distance inferred from dlon
        dy  : xarray.DataArray distance inferred from dlat
        """
        distance_1deg_equator = 2.0 * np.pi * radius * 1.0 / 360.0
        dxi = dlon * np.cos(lat / 180 * np.pi) * distance_1deg_equator
        dyi = ((lon * 0) + 1) * dlat * distance_1deg_equator
        return dxi, dyi

    if c.lon not in dsi.dims or c.lat not in dsi.dims:
        return dsi

    dsi = generate_grid_ds(dsi, {"X": c.lon, "Y": c.lat})
    xgcm_grid = xgcm.Grid(dsi, periodic=["X"])

    dlong = xgcm_grid.diff(dsi[c.lon], "X", boundary_discontinuity=360)
    dlonc = xgcm_grid.diff(dsi["lon_left"], "X", boundary_discontinuity=360)

    dlatg = xgcm_grid.diff(dsi[c.lat], "Y", boundary="fill", fill_value=np.nan)
    dlatc = xgcm_grid.diff(
        dsi["lat_left"], "Y", boundary="fill", fill_value=np.nan
    )

    dsi.coords["dxg"], dsi.coords["dyg"] = dll_dist(
        dlong, dlatg, dsi[c.lon], dsi[c.lat], radius
    )
    dsi.coords["dxc"], dsi.coords["dyc"] = dll_dist(
        dlonc, dlatc, dsi[c.lon], dsi[c.lat], radius
    )
    dsi.coords["area_c"] = dsi.dxc * dsi.dyc

    return dsi

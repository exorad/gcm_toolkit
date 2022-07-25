import os

import cubedsphere as cs
import cubedsphere.const as c
import numpy as np
import xgcm
from f90nml import Parser
from xgcm.autogenerate import generate_grid_ds


class MITgcmDataParser(Parser):
    def __init__(self):
        super().__init__()
        self.comment_tokens += '#'
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
    else:
        return default


def convert_winds_and_T(ds, T_dim, W_dim):
    """
    Convert winds and temperature in dataset.
    Winds are converted from Pa/s to m/s.
    Temperatures are converted from potential temperature to ordinary temperature

    Parameters
    ----------
    ds: Dataset
        dataset to be converted
    T_dim: str
        temperature datadimension to be converted
    W_dim: str
        vertical wind datadimension to be converted

    Returns
    -------
    ds: Dataset
        dataset with converted dimension
    """
    kappa = ds.attrs["R"] / ds.attrs["cp"]
    ds[T_dim] = ds[T_dim] * (ds[c.Z] / ds.attrs["p_ref"]) ** kappa

    # calculate scale height
    H = ds.attrs["R"] / ds.attrs["g"] * ds[T_dim]

    # calculate geometric height
    ds[c.Z_geo] = - H * np.log(ds[c.Z] / ds.attrs["p_ref"])

    if W_dim in ds:
        # interpolate vertical windspeed to cell center:
        if c.FACEDIM in ds.dims:
            grid = cs.init_grid_CS(ds=ds)
        else:
            grid = cs.init_grid_LL(ds=ds)

        W_interp = grid.interp(ds[W_dim], axis=c.Z, to="center")

        # convert vertical wind speed from Pa/s to m/s
        ds[W_dim] = - W_interp * H / ds[c.Z]

    return ds


def exorad_postprocessing(ds, outdir=None, datafile=None):
    """
    Preliminaray postprocessing on exorad dataset.
    This function converts the vertical windspeed from Pa into meters and saves attributes to the dataset.

    Parameters
    ----------
    ds: Dataset
        dataset to be extended
    outdir: string
        directory in which to find the data file (following the convention f'{outdir}/data')
    datafile: string
        alternatively specify datafile directly

    Returns
    ----------
    ds:
        Dataset to be returned
    """
    assert outdir is not None or datafile is not None, "please specify a datafile or a folder where we can find a datafile"

    if outdir is not None:
        datafile = f'{outdir}/data'

    # Add metadata
    radius = float(get_parameter(datafile, 'rSphere', 6370e3))  # planet radius in m
    P_rot = float(get_parameter(datafile, 'rotationperiod', 8.6164e4))  # planet rotationperiod in m
    attrs = {"p_ref": float(get_parameter(datafile, 'Ro_SeaLevel', 1.0e5)),  # bottom layer pressure in pascal
             "cp": float(get_parameter(datafile, 'atm_Cp', 1.004e3)),  # heat capacity at constant pressure
             "R": float(get_parameter(datafile, 'atm_Rd', 2.868571E2)),  # specific gas constant
             "g": float(get_parameter(datafile, 'gravity', 9.81)),  # surface gravity in m/s^2
             "dt": int(get_parameter(datafile, 'deltaT', 0.0)),  # time step size in s
             "R_p": radius,
             "P_rot": P_rot,
             "P_orb": P_rot,  # change, when this is available
             }

    ds.attrs.update(attrs)

    # Convert Temperature and winds
    if c.T in ds:
        ds = convert_winds_and_T(ds, c.T, c.W)
    if c.Ttave in ds:
        ds = convert_winds_and_T(ds, c.Ttave, c.wVeltave)

    # Add metrics to dataset
    if c.FACEDIM not in ds.dims:
        ds = add_distances(ds, radius=radius)

    return ds


def add_distances(ds, radius):
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
        dx = dlon * np.cos(lat/180*np.pi) * distance_1deg_equator
        dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
        return dx, dy

    if c.lon not in ds.dims or c.lat not in ds.dims:
        return ds

    ds = generate_grid_ds(ds, {'X': c.lon, 'Y': c.lat})
    xgcm_grid = xgcm.Grid(ds, periodic=['X'])

    dlong = xgcm_grid.diff(ds[c.lon], 'X', boundary_discontinuity=360)
    dlonc = xgcm_grid.diff(ds["lon_left"], 'X', boundary_discontinuity=360)

    dlatg = xgcm_grid.diff(ds[c.lat], 'Y', boundary='fill', fill_value=np.nan)
    dlatc = xgcm_grid.diff(ds["lat_left"], 'Y', boundary='fill', fill_value=np.nan)

    ds.coords['dxg'], ds.coords['dyg'] = dll_dist(dlong, dlatg, ds[c.lon], ds[c.lat], radius)
    ds.coords['dxc'], ds.coords['dyc'] = dll_dist(dlonc, dlatc, ds[c.lon], ds[c.lat], radius)
    ds.coords['area_c'] = ds.dxc * ds.dyc

    return ds

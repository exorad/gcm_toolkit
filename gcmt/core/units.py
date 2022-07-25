import astropy.units as u
import numpy as np

from .const import VARNAMES as c

ALLOWED_PUNITS = ['bar', 'Pa']
ALLOWED_TIMEUNITS = ['day', 'iter']


def convert_pressure(ds, current_unit, goal_unit):
    """
    Convert the pressure in dataset if current_unit and goal_unit differ.

    Parameters
    ----------
    ds: xarray.DataSet
        dataset for which the pressure should be converted
    current_unit: str
        Current unit
    goal_unit: str
        Unit to be changed to

    Returns
    -------
    ds: xarray.DataSet
        dataset with updated units
    """

    if current_unit not in ALLOWED_PUNITS or goal_unit not in ALLOWED_PUNITS:
        raise ValueError(f'current_unit and goal_unit need to be in {ALLOWED_PUNITS}')

    # Convert pressure
    for dim in {c['Z'], c['Z_l'], c['Z_p1'], c['Z_u']}:
        if dim in ds.dims:
            ds[dim] = (np.array(ds[dim]) * u.Unit(current_unit)).to(u.Unit(goal_unit)).value
    ds.attrs.update({'p_ref': (ds.p_ref * u.Unit(current_unit)).to(u.Unit(goal_unit)).value})

    # store the units in the dataset attributes
    ds.attrs['p_unit'] = goal_unit

    return ds


def convert_time(ds, current_unit, goal_unit):
    """
    Convert the time in dataset if current_unit and goal_unit differ.

    Parameters
    ----------
    ds: xarray.DataSet
        dataset for which the pressure should be converted
    current_unit: str
        Current unit
    goal_unit: str
        Unit to be changed to

    Returns
    -------
    ds: xarray.DataSet
        dataset with updated units
    """

    if current_unit not in ALLOWED_TIMEUNITS or goal_unit not in ALLOWED_TIMEUNITS:
        raise ValueError(f'current_unit and goal_unit need to be in {ALLOWED_TIMEUNITS}')

    if current_unit != goal_unit \
            and (current_unit == 'iter' or goal_unit == 'iter') \
            and ds.attrs.get('dt') is None:
        raise ValueError('You can only convert units if the timestep is known. ')

    # Convert time:
    if current_unit == 'iter' and goal_unit == 'day':
        ds[c["time"]] = ds.iter * ds.attrs["dt"] / (3600 * 24)
    if current_unit == 'day' and goal_unit == 'iter':
        ds[c["time"]] = ds.iter / ds.attrs["dt"] * (3600 * 24)

    # store the units in the dataset attributes
    ds.attrs['time_unit'] = goal_unit

    return ds

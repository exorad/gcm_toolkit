# ==============================================================
#                    GCM read in functionalities
# ==============================================================
#  This file contains all functionalities to read in from
#  different GCMs.
#
#  For MITgcm, we lean heavily on Aaron Schneider's cubedsphere
#  package, which in turn borrows functionality from xmitgcm and
#  xESMF and xgcm.
# ==============================================================

import glob
import os

import xarray as xr

import GCMtools.core.writer as wrt
from GCMtools.core.units import convert_pressure, convert_time


def m_read_from_mitgcm(gcmt, data_path, iters, d_lon=5, d_lat=4, loaded_ds=None, **kwargs):
    """
    Data read in for MITgcm output.

    Parameters
    ----------
    gcmt : GCMT
        GCMTools to which the data should be added
    data_path : str
        Folder path to the standard output of the GCM.
    iters : list, str
        The iteration (time step) of the input files to be read.
        If None, no data will be read.
        If 'last' (default), only the last iteration will be read.
        If 'all', all iterations will be read.
    data_file : str
        Full path to the 'data' input file of MITgcm. If None, the default
        location of the file is assumed to be: data_path/data
    **kwargs: dict
        pass down parameters to open_ascii_dataset from cubedsphere.

    Returns
    -------
    NoneType
        None
    """
    import cubedsphere as cs
    from .utils import exorad_postprocessing

    # determine the prefixes that should be loaded
    prefix = kwargs.pop('prefix', ["T", "U", "V", "W"])

    # determine the final iteration if needed
    if iters == 'last':
        all_iters = find_iters_mitgcm(data_path, prefix)
        iters = [max(all_iters)]
    elif iters == 'all':
        iters = find_iters_mitgcm(data_path, prefix)

    wrt.write_status('INFO', 'Iterations: ' + ", ".join([str(i) for i in iters]))

    if loaded_ds is not None:
        to_load = list(set(iters) - set(list(loaded_ds.iter.values)))
        if len(to_load) == 0:
            return loaded_ds
    else:
        to_load = iters

    # Currently, the read-in method is built using the wrapper functionality of
    # the cubedsphere package (Aaron Schneider)
    # see: https://cubedsphere.readthedocs.io/en/latest/index.html
    ds_ascii, grid = cs.open_ascii_dataset(data_path, iters=to_load, prefix=prefix, **kwargs)

    # regrid the dataset
    regrid = cs.Regridder(ds=ds_ascii, cs_grid=grid, d_lon=d_lon, d_lat=d_lat)
    ds = regrid()

    # convert wind, vertical dimension, time, ...
    ds = exorad_postprocessing(ds, outdir=data_path)

    convert_pressure(ds, current_unit='Pa', goal_unit=gcmt.p_unit)
    convert_time(ds, current_unit='iter', goal_unit=gcmt.time_unit)

    if loaded_ds is not None:
        ds = xr.merge([ds, loaded_ds])

    return ds


def find_iters_mitgcm(data_path, prefixes):
    """
    Helper method to list all iterations (time steps) that are present in the
    given MITgcm output directory.

    Parameters
    ----------
    data_path : str
        Folder path to the standard output of the GCM.

    Returns
    -------
    iterations : list of int
        List of all iterations that were found in the output folder.
    """
    iters_list = []
    for prefix in prefixes:
        files = glob.glob(os.path.join(data_path, "{}.*.data".format(prefix)))

        iters = [int(f.split('.')[-2]) for f in files]
        iters_list.append(iters)

    # find common iterations with data for all prefixes
    iterations = set.intersection(*[set(list) for list in iters_list])

    return iterations

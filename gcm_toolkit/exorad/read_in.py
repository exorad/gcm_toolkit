"""
==============================================================
                   GCM read in functionalities
==============================================================
 This file contains all functionalities to read in from
 different GCMs.

 For MITgcm, we lean heavily on Aaron Schneider's cubedsphere
 package, which in turn borrows functionality from xmitgcm and
 xESMF and xgcm.
==============================================================
"""
import glob
import os

import xarray as xr

from ..core import writer as wrt
from ..core.units import convert_pressure, convert_time


def m_read_from_mitgcm(
    tools,
    data_path,
    iters,
    exclude_iters=None,
    d_lon=5,
    d_lat=4,
    loaded_dsi=None,
    **kwargs,
):
    """
    Data read in for MITgcm output.

    Parameters
    ----------
    tools : GCMT
        GCMTools to which the data should be added
    data_path : str
        Folder path to the standard output of the GCM.
    iters : list, str
        The iteration (time step) of the input files to be read.
        If None, no data will be read.
        If 'last' (default), only the last iteration will be read.
        If 'all', all iterations will be read.
    exclude_iters: list, int, None
        List iterations that you don't want to load
        if None, no iterations from the list of iters will be excluded
        if list: exclude from iters
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
    prefix = kwargs.pop("prefix", ["T", "U", "V", "W"])

    # determine the final iteration if needed
    if iters == "last":
        all_iters = find_iters_mitgcm(data_path, prefix)
        iters = [max(all_iters)]
    elif iters == "all":
        iters = find_iters_mitgcm(data_path, prefix)

    if exclude_iters is not None:
        if isinstance(exclude_iters, int):
            exclude_iters = [exclude_iters]
        iters = set(iters) - set(exclude_iters)

    wrt.write_status(
        "INFO", "Iterations: " + ", ".join([str(i) for i in iters])
    )

    to_load = (
        list(set(iters) - set(list(loaded_dsi.iter.values)))
        if loaded_dsi is not None
        else iters
    )

    if len(to_load) == 0:
        if loaded_dsi is None:
            wrt.write_status("E-INFO", "No iterations selected to load.")
        return loaded_dsi

    # Currently, the read-in method is built using the wrapper functionality of
    # the cubedsphere package (Aaron Schneider)
    # see: https://cubedsphere.readthedocs.io/en/latest/index.html
    dsi_ascii, grid = cs.open_ascii_dataset(
        data_path, iters=to_load, prefix=prefix, **kwargs
    )

    # regrid the dataset
    filename = (  # generate random filename for weights to be deleted
        "tmp_gcmt_reg_weights_xya"
    )
    regrid = cs.Regridder(
        ds=dsi_ascii, cs_grid=grid, d_lon=d_lon, d_lat=d_lat, filename=filename
    )
    _ = [
        os.remove(f) for f in glob.glob(filename + "*.nc")
    ]  # delete aux weights
    dsi = regrid()

    # convert wind, vertical dimension, time, ...
    dsi = exorad_postprocessing(dsi, outdir=data_path)

    convert_pressure(dsi, current_unit="Pa", goal_unit=tools.p_unit)
    convert_time(dsi, current_unit="iter", goal_unit=tools.time_unit)

    if loaded_dsi is not None:
        dsi = xr.merge([dsi, loaded_dsi])

    return dsi


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
        files = glob.glob(os.path.join(data_path, f"{prefix}.*.data"))

        iters = [int(f.split(".")[-2]) for f in files]
        iters_list.append(iters)

    # find common iterations with data for all prefixes
    iterations = set.intersection(*[set(list) for list in iters_list])

    return iterations

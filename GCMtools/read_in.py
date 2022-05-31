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

import os
import cubedsphere as cs

def m_read_from_mitgcm(gcmt, data_path, iters):
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

    Returns
    -------
    NoneType
        None
    """
    # determine the final iteration if needed
    if iters == 'last':
        all_iters = find_iters_mitgcm(data_path)
        iters = [max(all_iters)]

    print('[INFO] -- Preparing to read from MITgcm data directory:' + data_path)
    print('          Iterations: ' + ", ".join([str(i) for i in iters]))

    # Currently, the read-in method is built using the wrapper functionality of
    # the cubedsphere package (Aaron Schneider)
    # see: https://cubedsphere.readthedocs.io/en/latest/index.html
    ds_ascii, grid = cs.open_ascii_dataset(data_path, iters=iters, prefix = ["T","U","V","W"])

    # regrid the dataset
    regrid = cs.Regridder(ds_ascii, grid)
    ds = regrid()

    # convert wind, vertical dimension, time, ...
    ds = cs.exorad_postprocessing(ds, outdir=data_path)

    return ds


def find_iters_mitgcm(data_path):
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
    iterations = []

    for data_file in os.listdir(data_path):
        # look for files matching the pattern T.xxxxxxxxxx.data
        if data_file[:2] == 'T.' and data_file[-5:] == '.data':
            # if it matches, add the iteration number to the list
            iterations.append(int(data_file[2:-5]))

    return iterations

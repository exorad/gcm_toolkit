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

def m_read_from_mitgcm(gcmt, data_path, iters, data_file=None):
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

    Returns
    -------
    NoneType
        None
    """
    # determine the final iteration if needed
    if iters == 'last':
        all_iters = find_iters_mitgcm(data_path)
        iters = [max(all_iters)]

    print('[INFO] Preparing to read from MITgcm data directory:' + data_path)
    print('       Iterations: ' + ", ".join([str(i) for i in iters]))

    # Currently, the read-in method is built using the wrapper functionality of
    # the cubedsphere package (Aaron Schneider)
    # see: https://cubedsphere.readthedocs.io/en/latest/index.html
    ds_ascii, grid = cs.open_ascii_dataset(data_path, iters=iters, prefix = ["T","U","V","W"])

    # regrid the dataset
    regrid = cs.Regridder(ds_ascii, grid)
    ds = regrid()

    # convert wind, vertical dimension, time, ...
    ds = cs.exorad_postprocessing(ds, outdir=data_path)

    # set location of data-file
    if data_file is None:
        data_file = data_path+'/data'

    # supplement the dataset attributes with necessary values
    ds.attrs['R_p'] = ds.attrs.pop('radius') # rename radius
    P_rot = float(get_data_parameter(data_file, 'rotationperiod')) # rotation period in seconds
    ds.attrs['P_rot'] = P_rot
    ds.attrs['P_orb'] = P_rot # TODO : update when non-synchronously rotating planets are allowed

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

    for f in os.listdir(data_path):
        # look for files matching the pattern T.xxxxxxxxxx.data
        if f[:2] == 'T.' and f[-5:] == '.data':
            # if it matches, add the iteration number to the list
            iterations.append(int(f[2:-5]))

    return iterations

def get_data_parameter(data_file, keyword):
    """
    Function to parse the MITgcm 'data' file and return the parameter values
    of the given specific keyword.

    Parameters
    ----------
    data_file: string
        Full path to the MITgcm data file.
    keyword: string
        Parameter of which the value is required.

    Returns
    -------
    value: string
        The value associated with the given keyword is returned as a string (!).
    """
    # Check if the data file exists
    if not os.path.isfile(data_file):
        raise FileNotFoundError('[ERROR] Data-file not found at: '+data_file)
    else: # if it does:
        valueFound = False
        with open(data_file, 'r') as f:
            for line in f.readlines():
                if keyword in line and line[0] != '#': # find uncommented line with the keyword
                    line = line.replace(" ", "")       # remove all unnecessary spaces
                    value = line.split(keyword+'=', 1)[1].rstrip().rstrip(',') # separate keyword and value
                    valueFound = True
                    return value
        if not valueFound:
            raise NameError(r'[ERROR] No value could be associated with the keyword: '+keyword)

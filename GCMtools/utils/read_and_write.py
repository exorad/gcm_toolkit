import glob
import os

import numpy as np
import xarray as xr

import GCMtools.core.writer as wrt
from GCMtools.core.const import VARNAMES as c
from GCMtools.core.units import convert_time, convert_pressure
from GCMtools.utils.passport import is_the_data_basic


def m_read_raw(gcmt, gcm, data_path, iters='last', load_existing=False, tag=None, **kwargs):
    """
    General read in function for GCM data

    Parameters
    ----------
    gcm : str
        Type of GCM, must be 'MITgcm'.
    data_path : str
        Folder path to the standard output of the GCM.
    iters : list, str
        The iteration (time step) of the input files to be read.
        If None, no data will be read.
        If 'last' (default), only the last iteration will be read.
        If 'all', all iterations will be read.
    load_existing: bool
        Set to false if you want to overwrite already loaded data
        Set to true if you want to increment already loaded data
    tag : str
        Tag to reference the simulation in the collection of models.
    kwargs: dict
        Additional options passed down to read functions
    """

    # call the required GCM read-in method
    if gcm == 'MITgcm':
        wrt.write_status('STAT', 'Read in raw MITgcm data')
        wrt.write_status('INFO', 'File path: ' + data_path)
        from GCMtools.exorad import m_read_from_mitgcm

        if tag is not None and load_existing:
            try:
                loaded_ds = gcmt.get_models(tag)
            except KeyError:
                loaded_ds = None
        else:
            loaded_ds = None

        ds = m_read_from_mitgcm(gcmt, data_path, iters, loaded_ds=loaded_ds, **kwargs)
    else:
        wrt.write_status('ERROR', 'The selected GCM type "' + gcm + '" is not supported')

    _add_attrs_and_store(gcmt, ds, tag)

def m_read_reduced(gcmt, data_path, tag=None, time_unit_in='iter', p_unit_in='Pa'):
    """
    Read in function for GCM data that has been reduced and saved according
    to the GCMtools GCMDataset format.

    Parameters
    ----------
    data_path : str
        Folder path to the reduced (GCMtools) data.
    time_unit_in: str
        units of time dimension in input dataset
    p_unit_in: str
        units of pressure dimensions in input dataset
    tag : str
        Tag to reference the simulation in the collection of models.
    """
    # print information
    wrt.write_status('STAT', 'Read in reduced data')
    wrt.write_status('INFO', 'File path: ' + data_path)
    wrt.write_status('INFO', 'Time unit of data: ' + time_unit_in)
    wrt.write_status('INFO', 'pressure unit of data: ' + p_unit_in)

    # read dataset using xarray functionalities
    ds = xr.open_dataset(data_path)

    ds = convert_time(ds, current_unit=time_unit_in, goal_unit=gcmt.time_unit)
    ds = convert_pressure(ds, current_unit=p_unit_in, goal_unit=gcmt.p_unit)

    _add_attrs_and_store(gcmt, ds, tag)

    wrt.write_status('INFO', 'Tag: ' + tag)


def _add_attrs_and_store(gcmt, ds, tag):
    # if no tag is given, models are just numbered as they get added
    if tag is None:
        tag = str(len(gcmt.get_models(always_dict=True)))
        print('[WARN] -- No tag provided. This model is stored with tag: ' + tag)

    wrt.write_status('INFO', 'Tag: ' + tag)

    # store tag in the dataset attributes
    ds.attrs['tag'] = tag
    # check if the dataset has all necessary GCMtools attributes
    if not is_the_data_basic(ds):
        raise ValueError('This dataset is not supported by GCMtools\n')

    # store dataset
    gcmt[tag] = ds


def m_save(gcmt, dir, method='nc', update_along_time=False, tag=None):
    """
    Save function to store current member variables.

    Parameters
    ----------
    dir : str
        directory at which the GCMtools datasets should be stored.
    method : str, optional
        Datasets can be stored as '.zarr' or '.nc'. Decide which type you prefer.
        Defaults to '.nc'.
    update_along_time: str, optional
        Decide if you want to update already saved datasets along the timedimension.
        This only works with method='zarr'.
    tag: str, optional
        tag of the model that should be loaded. Will save all available models by default.

    Returns
    -------
    NoneType
        None
    """

    # print information
    wrt.write_status('STAT', 'Save current GCMs within GCMtools')
    wrt.write_status('INFO', 'File path: ' + dir)
    if tag is None:
        wrt.write_message('INFO', 'Tag: All tags were stored')
    else:
        wrt.write_status('INFO', 'Tag: ' + tag)
    wrt.write_status('INFO', 'method: ' + method)
    wrt.write_status('INFO', 'Update old data?: ' + str(update_along_time))

    if method not in ['nc', 'zarr']:
        raise NotImplementedError("Please use zarr or nc.")

    for key, model in gcmt.get_models(always_dict=True).items():
        if tag is not None and tag != key:
            continue

        filename = os.path.join(dir, f"{key}.{method}")

        if method == 'nc':
            if os.path.isfile(filename):
                os.remove(filename)
            model.to_netcdf(filename)
        elif method == 'zarr':
            if os.path.isdir(filename) and update_along_time:
                # Relies on https://stackoverflow.com/questions/65339851/xarray-dataset-to-zarr-overwrite-data-if-exists-with-append-dim
                # read structure of dataset to see what's on disk
                ds_ondisk = xr.open_zarr(filename)

                # get index of first new datapoint
                start_ix, = np.nonzero(~np.isin(model[c['time']], ds_ondisk[c['time']]))

                if len(start_ix) > 0:
                    # region of new data
                    region_new = slice(start_ix[0], model[c['time']].size)

                    # append structure of new data (compute=False means no data is written)
                    model.isel(time=region_new).to_zarr(filename, append_dim=c['time'], compute=True)
            else:
                model.to_zarr(filename, mode='w')


def m_load(gcmt, dir, method='nc', tag=None):
    """
    Load function to load stored member variables.

    Parameters
    ----------
    dir : str
        directory at which the GCMtools datasets are stored
    method : str, optional
        Should be the same method with which you stored the data
    tag: str, optional
        tag of the model that should be loaded. Will load all available models by default.
    """

    # print information
    wrt.write_status('STAT', 'Load saved GCMs to GCMtools')
    wrt.write_status('INFO', 'File path: ' + dir)
    if tag is None:
        wrt.write_message('INFO', 'Tag: All tags were stored')
    else:
        wrt.write_status('INFO', 'Tag: ' + tag)
    wrt.write_status('INFO', 'method: ' + method)

    if method not in ['nc', 'zarr']:
        raise NotImplementedError("Please use zarr or nc.")

    if tag is None:
        available_datasets = glob.glob(f'{dir}/*.{method}')
    else:
        available_datasets = glob.glob(f'{dir}/{tag}.{method}')

    if len(available_datasets) == 0:
        print(f'[INFO] No data available to load for method {method}')

    for file in available_datasets:
        head, tail = os.path.split(file)
        tag = tail.replace(f'.{method}', '')
        if method == 'zarr':
            ds = xr.open_zarr(file)
        elif method == 'nc':
            ds = xr.open_dataset(file)

        ds = convert_time(ds, current_unit=ds.attrs.get('time_unit'), goal_unit=gcmt.time_unit)
        ds = convert_pressure(ds, current_unit=ds.attrs.get('p_unit'), goal_unit=gcmt.p_unit)

        gcmt[tag] = ds

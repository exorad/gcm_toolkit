# ==============================================================
#                       GCMtools Main Class
# ==============================================================
#  This class is the user interace of the GCMtools
#  functionalities. The goal is to have a clean and easy to use
#  environment for new users of GCMs while allowing direct
#  access to the data for more experienced users.
# ==============================================================
import os.path
import glob

import xarray as xr
import numpy as np
from collections import UserDict

from .passport import is_the_data_basic


class GCMT:
    """
    The main GCMtools class with which the user can interact.

    Attributes
    ----------
    models : dict
        Dictionary containing all of the 3D GCM models that are stored in the
        memory, with their respective tags.

    Methods
    -------
    read_raw(self, gcm, data_path, tag=None):
        Read in the raw data output from a GCM.
    read_reduced(self, data_path, tag=None):
        Read in the previously reduced GCM
    """

    def __init__(self, use_bar=True, use_days=True):
        """
        Constructor for the GCMtools class.

        Parameters
        ----------
        models : dict
            Dictionary containing all of the GCM datasets that have been read.
        """

        # Initialize empty dictionary to store all GCM models
        self._models = GCMDatasetCollection()
        self.use_bar = use_bar
        self.use_days = use_days

    @property
    def models(self):
        if len(self._models) > 1:
            return self._models
        else:
            return list(self._models.values())[0]

    def get_models(self, tag=None):
        """
        Function return all GCMs in memory. If a tag is given, only return this
        one.

        Parameters
        ----------
        tag : str
            Name of the model that should be returned.

        Returns
        -------
        selected_models : GCMDatasetCollection
            All models in self._models, or only the one with the right tag.
        """
        # If no tag is given, return all models
        if tag is None:
            return self.models
        # If the tag is valid, return the corresponding model
        elif isinstance(tag, str):
            return self._models[tag]
        # If the tag is not a string, raise an error
        else:
            raise ValueError('The given tag is not a string.')

    def read_raw(self, gcm, data_path, iters='last', load_existing=False, tag=None, **kwargs):
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
            from GCMtools.exorad import m_read_from_mitgcm

            if tag is not None and load_existing:
                loaded_ds = self._models.get(tag)
            else:
                loaded_ds = None

            ds = m_read_from_mitgcm(self, data_path, iters, loaded_ds=loaded_ds, **kwargs)
        else:
            raise ValueError('The selected GCM type "' + gcm +
                             '" is not supported')

        self._add_attrs_and_store(ds, tag)

    def read_reduced(self, data_path, tag=None):
        """
        Read in function for GCM data that has been reduced and saved according
        to the GCMtools GCMDataset format.

        Parameters
        ----------
        data_path : str
            Folder path to the reduced (GCMtools) data.
        tag : str
            Tag to reference the simulation in the collection of models.
        """
        # read dataset using xarray functionalities
        ds = xr.open_dataset(data_path)

        self._add_attrs_and_store(ds, tag)

    def _add_attrs_and_store(self, ds, tag):
        # if no tag is given, models are just numbered as they get added
        if tag is None:
            print('[WARN] -- No tag provided. This model is stored with tag: ' + str(len(self._models)))
            tag = str(len(self._models))
        # store tag in the dataset attributes
        ds.attrs['tag'] = tag
        # store the units in the dataset attributes
        ds.attrs['p_unit'] = 'bar' if self.use_bar else 'Pa'
        ds.attrs['time_unit'] = 'days' if self.use_days else 'iter'
        # check if the dataset has all necessary GCMtools attributes
        if not is_the_data_basic(ds):
            raise ValueError('This dataset is not supported by GCMtools\n')

        # store dataset
        self._models[tag] = ds

    def save(self, dir, update_along_time = False, tag = None):
        """
        Save function to store current member variables.

        Relies on https://stackoverflow.com/questions/65339851/xarray-dataset-to-zarr-overwrite-data-if-exists-with-append-dim

        Parameters
        ----------
        example : int
            Short description of the variable.
        update_along_time: bool
            Append already written files along time dimension (less IO)

        Returns
        -------
        NoneType
            None
        """
        for key, model in self._models.items():
            if tag is not None and tag != key:
                continue

            filename = os.path.join(dir, f"{key}.zarr")

            if os.path.isdir(filename) and update_along_time:
                # read structure of dataset to see what's on disk
                ds_ondisk = xr.open_zarr(filename)

                # get index of first new datapoint
                start_ix, = np.nonzero(~np.isin(model.time, ds_ondisk.time))

                if len(start_ix)>0:
                    # region of new data
                    region_new = slice(start_ix[0], model.time.size)

                    # append structure of new data (compute=False means no data is written)
                    model.isel(time=region_new).to_zarr(filename, append_dim='time', compute=True)
            else:
                model.to_zarr(filename, mode='w')


    def load(self, dir, tag = None):
        """
        Save function to store current member variables.

        Parameters
        ----------
        example : int
            Short description of the variable.

        Returns
        -------
        NoneType
            None
        """
        if tag is None:
            available_datasets = glob.glob(f'{dir}/*.zarr')
        else:
            available_datasets = glob.glob(f'{dir}/{tag}.zarr')

        for file in available_datasets:
            head, tail = os.path.split(file)
            tag = tail.replace('.zarr', '')
            self._models[tag] = xr.open_zarr(file)


# ------------------------------------------------------------------------------


class GCMDatasetCollection(UserDict):
    """
    This class represents a collection of 3D GCM Datasets.
    A GCMDatasetCollection is a dictionary with in which GCM models are loaded
    with a tag.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # TODO add additional functionalitiess

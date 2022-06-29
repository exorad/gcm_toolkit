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
from .units import ALLOWED_PUNITS, ALLOWED_TIMEUNITS, convert_time, convert_pressure
from .const import VARNAMES as c
import GCMtools.gcm_plotting as gcmplt


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

    def __init__(self, p_unit='bar', time_unit='day'):
        """
        Constructor for the GCMtools class.

        Parameters
        ----------
        p_unit: str, optional
            Set the unit that is used internally for pressure related things
        time_unit: str, optional
            Set the unit that is used internally for time related things
        """

        # Initialize empty dictionary to store all GCM models
        self._models = GCMDatasetCollection()

        if p_unit not in ALLOWED_PUNITS:
            raise ValueError(f"Please use a pressure unit from {ALLOWED_PUNITS}")
        self.p_unit = p_unit

        if time_unit not in ALLOWED_TIMEUNITS:
            raise ValueError(f"Please use a time unit from {ALLOWED_TIMEUNITS}")
        self.time_unit = time_unit

    @property
    def models(self):
        """
        returns dataset if only one model is loaded, otherwise collection of datasets.
        """
        if len(self._models) > 1:
            return self._models
        else:
            return list(self._models.values())[0]

    def _get_one_model(self, tag=None):
        """
        Helper Function that raises an error, if more than one model is selected.

        Parameters
        ----------
        tag: str, optional
            Name of the model that should be returned
        Returns
        -------
        ds: xarray Dataset
            Selected model
        """

        # select the appropriate dataset
        ds = self.get_models(tag=tag)
        # if a collection is given (because multiple datasets are available, and
        # the tag is not provided), avoid ambiguity by raising an error
        if isinstance(ds, GCMDatasetCollection) and len(ds) > 1 and tag is None:
            raise RuntimeError('Ambiguous task. Please provide a tag.')

        return ds

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
        selected_models : GCMDatasetCollection or xarray Dataset
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

    def read_reduced(self, data_path, tag=None, time_unit_in='iter', p_unit_in='Pa'):
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
        # read dataset using xarray functionalities
        ds = xr.open_dataset(data_path)

        ds = convert_time(ds, current_unit=time_unit_in, goal_unit=self.time_unit)
        ds = convert_pressure(ds, current_unit=p_unit_in, goal_unit=self.p_unit)

        self._add_attrs_and_store(ds, tag)

    def _add_attrs_and_store(self, ds, tag):
        # if no tag is given, models are just numbered as they get added
        if tag is None:
            print('[WARN] -- No tag provided. This model is stored with tag: ' + str(len(self._models)))
            tag = str(len(self._models))
        # store tag in the dataset attributes
        ds.attrs['tag'] = tag
        # check if the dataset has all necessary GCMtools attributes
        if not is_the_data_basic(ds):
            raise ValueError('This dataset is not supported by GCMtools\n')

        # store dataset
        self._models[tag] = ds

    def save(self, dir, method='nc', update_along_time=False, tag=None):
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
        if method not in ['nc', 'zarr']:
            raise NotImplementedError("Please use zarr or nc.")

        for key, model in self._models.items():
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

    def load(self, dir, method='nc', tag=None):
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

            ds = convert_time(ds, current_unit=ds.attrs.get('time_unit'), goal_unit=self.time_unit)
            ds = convert_pressure(ds, current_unit=ds.attrs.get('p_unit'), goal_unit=self.p_unit)

            self._models[tag] = ds

    def add_horizontal_average(self, var_key, var_key_out=None, area_key='area_c', tag=None):
        """
        Calculate horizontal averaged quantities. Horizontal averages
        are calculated as:

            \bar q = \int{q dA}/\int{dA}

        Parameters
        ----------
        var_key: str
            The key of the variable quantity that should be averaged.
        tag : str, optional
            The tag of the dataset that should be used. If no tag is provided,
            and multiple datasets are available, an error is raised.
        var_key_out: str, optional
            variable name used to store the outcome. If not provided, this script will just
            return the averages and not change the dataset inplace.
        area_key: str, optional
            Variable key in the dataset for the area of grid cells
        """
        ds = self._get_one_model(tag)
        avg = (ds[area_key]*ds[var_key]).sum(dim=[c['lon'],c['lat']])/ds[area_key].sum(dim=[c['lon'],c['lat']])

        if var_key_out is not None:
            ds.update({var_key_out: avg})

        return avg

    def add_meridional_overturning(self, v_data='V', var_key_out=None, tag=None):
        """
        Calculate meridional overturning streamfunction.

        See e.g. Carone et al. (2017), Eq. 7

        Parameters
        ----------
        var_key: str
            The key of the meridional velocity that should be used to calculate the outcome
        tag : str, optional
            The tag of the dataset that should be used. If no tag is provided,
            and multiple datasets are available, an error is raised.
        var_key_out: str, optional
            variable name used to store the outcome. If not provided, this script will just
            return the overturning circulation and not change the dataset inplace.
        """
        ds = self._get_one_model(tag)
        V_integral = ds[v_data].cumulative_integrate(coord='Z')

        if ds.attrs.get('p_unit') == 'bar':
            # convert to SI, if needed
            V_integral = V_integral / 1.0e5

        psi = 2*np.pi*np.cos(ds.lat/180*np.pi)*ds.R_p/ds.g*V_integral

        if var_key_out is not None:
            ds.update({var_key_out: psi})

        return psi

    def isobaric_slice(self, var_key, p, tag=None, **kwargs):
        """
        Plot an isobaric slice of the given quantity at the given pressure
        level. The user can specify the DataSet to be plotted by providing the
        corresponding tag.

        Parameters
        ----------
        var_key : str
            The key of the variable quantity that should be plotted.
        p : float
            Pressure level for the isobaric slice to be plotted, expressed in
            the units specified in the dataset attributes (e.g., init of GCMT object).
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        # select the appropriate dataset
        ds = self._get_one_model(tag)
        return gcmplt.isobaric_slice(ds, var_key, p, **kwargs)

    def time_evol(self, var_key, tag=None, **kwargs):
        """
        Function that plots the time evolution of a quantity in a 1D line collection plot, where the colorscale can be related to the time evolution.
        Note: var_key needs to contain data that is 2D in time and pressure.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        ds = self._get_one_model(tag)
        return gcmplt.time_evol(ds, var_key, **kwargs)

    def zonal_mean(self, var_key, tag=None, **kwargs):
        """
        Plot a zonal mean average of a quantity for the given dataset.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        # select the appropriate dataset
        ds = self._get_one_model(tag)
        return gcmplt.zonal_mean(ds, var_key, **kwargs)


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

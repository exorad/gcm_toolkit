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

# import GCMDatasetCollection class
from GCMtools.GCMDatasetCollection import GCMDatasetCollection

# import core functionalities
from GCMtools.core.units import ALLOWED_PUNITS, ALLOWED_TIMEUNITS, convert_time, convert_pressure
from GCMtools.core.const import VARNAMES as c

# import utility functions
import GCMtools.utils.read_and_write as raw
import GCMtools.utils.gcm_plotting as gcmplt
import GCMtools.utils.manipulations as mani


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

    # =============================================================================================================
    #   Data handling
    # =============================================================================================================

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
        return self._models.get_models()

    # =============================================================================================================
    #   Data manipulation
    # =============================================================================================================
    def add_horizontal_average(self, var_key, var_key_out=None, area_key='area_c', tag=None):
        """
        Calculate horizontal averaged quantities. Horizontal averages
        are calculated as:

            \bar q = \int{q dA}/\int{dA}

        Parameters
        ----------
        var_key: str
            The key of the variable quantity that should be plotted.
        tag : str, optional
            The tag of the dataset that should be used. If no tag is provided,
            and multiple datasets are available, an error is raised.
        var_key_out: str, optional
            variable name used to store the outcome. If not provided, this script will just
            return the averages and not change the dataset inplace.
        area_key: str, optional
            Variable key in the dataset for the area of grid cells

        Returns
        -------
        TODO
        """
        return mani.m_add_horizontal_average(self, var_key, var_key_out=None, area_key='area_c', tag=None)

    # =============================================================================================================
    #   Reading and writing functions
    # =============================================================================================================
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
        return raw.m_read_raw(self, gcm, data_path, iters=iters, load_existing=load_existing, tag=tag, **kwargs)

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
        return raw.m_read_reduced(self, data_path, tag=None, time_unit_in='iter', p_unit_in='Pa')

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
        return raw.m_save(self, dir, method='nc', update_along_time=False, tag=None)

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
        return raw.m_load(self, dir, method='nc', tag=None)

    # =============================================================================================================
    #   Plotting functions
    # =============================================================================================================

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
        ds = self._models.get_one_model(tag)
        return gcmplt.isobaric_slice(ds, var_key, p, **kwargs)

    def time_evol(self, var_key, tag=None, **kwargs):
        """
        Function that plots the time evolution of a quantity in a 1D line collection plot,
        where the colorscale can be related to the time evolution.
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
        ds = self._models.get_one_model(tag)
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
        ds = self._models.get_one_model(tag)
        return gcmplt.zonal_mean(ds, var_key, **kwargs)

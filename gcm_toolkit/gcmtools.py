"""
==============================================================
                      gcm_toolkit Main Class
==============================================================
 This class is the user interace of the gcm_toolkit
 functionalities. The goal is to have a clean and easy to use
 environment for new users of GCMs while allowing direct
 access to the data for more experienced users.
==============================================================
"""
import xarray

from .core import writer as wrt
from .utils import gcm_plotting as gcmplt
from .utils import manipulations as mani
from .utils import read_and_write as raw
from .gcm_dataset_collection import GCMDatasetCollection
from .core.units import ALLOWED_PUNITS, ALLOWED_TIMEUNITS
from .utils.passport import is_the_data_basic
from .utils.interface import PrtInterface


class GCMT:
    """
    The main gcm_toolkit class with which the user can interact.

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

    def __init__(self, p_unit="bar", time_unit="day", write="on"):
        """
        Constructor for the gcm_toolkit class.

        Parameters
        ----------
        p_unit: str, optional
            Set the unit that is used internally for pressure related things
        time_unit: str, optional
            Set the unit that is used internally for time related things
        write: str, optional
            Set the level of logging information. Check core.writer for more infos
        """

        # Initialize empty dictionary to store all GCM models
        self._models = GCMDatasetCollection()

        # check units
        if p_unit not in ALLOWED_PUNITS:
            raise ValueError(
                f"Please use a pressure unit from {ALLOWED_PUNITS}"
            )
        self.p_unit = p_unit

        if time_unit not in ALLOWED_TIMEUNITS:
            raise ValueError(
                f"Please use a time unit from {ALLOWED_TIMEUNITS}"
            )
        self.time_unit = time_unit

        # initialize writing function (to file or to console)
        wrt.writer_setup(write)

        # print welcome message and information
        wrt.write_hline()
        hello = "Welcome to gcm_toolkit"
        wrt.write_message(
            hello,
            color="WARN",
            spacing=(wrt.Writer.line_length - len(hello)) // 2,
        )
        wrt.write_hline()
        wrt.write_status("STAT", "Set up gcm_toolkit")
        wrt.write_status("INFO", "pressure units: " + self.p_unit)
        wrt.write_status("INFO", "time units: " + self.time_unit)

    # ==============================================================================================
    #   Data handling
    # ==============================================================================================

    @property
    def models(self):
        """
        Shortcut to return all GCMs in memory.

        Returns
        -------
        selected_models : GCMDatasetCollection or xarray Dataset
            All models in self._models
        """
        return self.get_models()

    def __iter__(self):
        """
        Returns an iterator of GCMT
        """
        return iter(self.get_models(always_dict=True).items())

    def __getitem__(self, tag):
        """
        Access models of this dataset

        Parameters:
        -----------
        tag: str
            Tag of the model that should be returned

        Returns
        -------
        selected_model : xarray Dataset
        """
        if tag not in self._models.keys():
            raise ValueError(
                "The provided tag does not exist in the collection"
            )
        if isinstance(tag, str):
            return self.get_one_model(tag=tag, raise_error=True)

    def __setitem__(self, tag, dsi):
        """
        Add or replaces a dataset

        Parameters
        ----------
        tag: str
           Tag at which the model should be stored
        dsi: xarray Dataset
           Dataset to add
        """
        self._replace_model(tag, dsi)

    def __len__(self) -> int:
        return len(self._models)

    def __bool__(self) -> bool:
        return bool(self._models)

    def get(self, tag, default=None):
        """
        Pythonic get, will return default if tag is not available.

        Parameters
        ----------
        tag: str
            Name of the model that should be returned
        default: Any
            Default value in case that the tag is not available in the dataset

        Returns
        -------
        dsi or default, depending on whether dsi is available
        """
        if dsi := self.get_one_model(tag, raise_error=False):
            return dsi
        return default

    def get_one_model(self, tag=None, raise_error=True):
        """
        Like get_models, but raises an error or return None,
        if more than one model is selected.

        Parameters
        ----------
        tag: str, optional
            Name of the model that should be returned
        raise_error: bool, optional
            If true, function will raise error, else will return
            None if not one model is selected

        Returns
        -------
        dsi: xarray Dataset
            Selected model
        """
        return self._models.get_one_model(tag, raise_error)

    def get_models(self, tag=None, always_dict=False):
        """
        Function return all GCMs in memory. If a tag is given, only return this
        one.

        Parameters
        ----------
        tag : str
            Name of the model that should be returned.
        always_dict: bool
            Force result to be a dictionary (if tag is None)

        Returns
        -------
        selected_models : GCMDatasetCollection or xarray Dataset
            All models in self._models, or only the one with the right tag.
            Will definetly be GCMDatasetCollection if always_dict=True
        """
        return self._models.get_models(tag, always_dict)

    def _replace_model(self, tag, dsi):
        """
        Add or replaces a dataset. Do some checks beforehand.

        Parameters
        ----------
        tag: str
           Tag at which the model should be stored
        dsi: xarray Dataset
           Dataset to add
        """
        if not isinstance(tag, str):
            raise ValueError("The provided tag needs to be a string.")
        if not isinstance(dsi, xarray.Dataset):
            raise ValueError(
                "The provided input dataset needs to be a xarray Dataset."
            )
        if not is_the_data_basic(dsi):
            raise ValueError("The provided input dataset is not compatible.")

        if dsi.attrs.get("p_unit") != self.p_unit:
            raise NotImplementedError(
                "Unit conversion is needed and not yet implemented"
            )
        if dsi.attrs.get("time_unit") != self.time_unit:
            raise NotImplementedError(
                "Unit conversion is needed and not yet implemented"
            )

        dsi.attrs.update({"tag": tag})
        self._models[tag] = dsi

    # ==============================================================================================
    #   Data manipulation
    # ==============================================================================================
    def add_horizontal_average(
        self, var_key, var_key_out=None, area_key="area_c", tag=None
    ):
        """
        Calculate horizontal averaged quantities. Horizontal averages
        are calculated as area-weighted quantities q, so that:

            bar{q} = int{q dA}/int{dA}

        Parameters
        ----------
        var_key: str
            The key of the variable quantity that should be plotted.
        var_key_out: str, optional
            variable name used to store the outcome.
            If not provided, this script will just
            return the averages and not change the dataset inplace.
        area_key: str, optional
            Variable key in the dataset for the area of grid cells
        tag : str, optional
            The tag of the dataset that should be used.
            If no tag is provided, and multiple datasets are available,
            an error is raised.

        Returns
        -------
        avg: xarray.DataArray
            Averaged data

        """
        dsi = self.get_one_model(tag)
        return mani.m_add_horizontal_average(
            dsi, var_key, var_key_out=var_key_out, area_key=area_key
        )

    def add_meridional_overturning(
        self, v_data="V", var_key_out=None, tag=None
    ):
        """
        Calculate meridional overturning streamfunction.

        See e.g. Carone et al. (2017), Eq. 7

        Parameters
        ----------
        v_data: str
            The key that holds the data to meridional velocity
        tag : str, optional
            The tag of the dataset that should be used.
             If no tag is provided,
            and multiple datasets are available, an error is raised.
        var_key_out: str, optional
            variable name used to store the outcome.
            If not provided, this script will just
            return the overturning circulation
            and not change the dataset inplace.

        Returns
        -------
        psi: xarray.DataArray
            Horizontal overturning
        """
        dsi = self.get_one_model(tag)
        return mani.m_add_meridional_overturning(
            dsi, v_data=v_data, var_key_out=var_key_out
        )

    # ======================================================
    #   Reading and writing functions
    # ======================================================
    def read_raw(
        self,
        gcm,
        data_path,
        iters="last",
        load_existing=False,
        tag=None,
        **kwargs,
    ):
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
        return raw.m_read_raw(
            self,
            gcm,
            data_path,
            iters=iters,
            load_existing=load_existing,
            tag=tag,
            **kwargs,
        )

    def read_reduced(
        self, data_path, tag=None, time_unit_in=None, p_unit_in=None
    ):
        """
        Read in function for GCM data that has been reduced a
        nd saved according to the gcm_toolkit GCMDataset format.

        Parameters
        ----------
        data_path : str
            Folder path to the reduced (gcm_toolkit) data.
        time_unit_in: str, None
            units of time dimension in input dataset.
            If None, try to read from nc file (ds.attr.time_unit)
        p_unit_in: str, None
            units of pressure dimensions in input dataset,
            If None, try to read from nc file (ds.attr.p_unit)
        tag : str
            Tag to reference the simulation in the collection of models.
        """
        return raw.m_read_reduced(
            self,
            data_path,
            tag=tag,
            time_unit_in=time_unit_in,
            p_unit_in=p_unit_in,
        )

    def save(self, direct, method="nc", update_along_time=False, tag=None):
        """
        Save function to store current member variables.

        Parameters
        ----------
        direct : str
            directory at which the gcm_toolkit datasets should be stored.
        method : str, optional
            Datasets can be stored as '.zarr' or '.nc'.
            Decide which type you prefer. Defaults to '.nc'.
        update_along_time: bool, optional
            Decide if you want to update already saved datasets
            along the timedimension. This only works with method='zarr'.
        tag: str, optional
            tag of the model that should be loaded.
            Will save all available models by default.

        Returns
        -------
        NoneType
            None
        """
        return raw.m_save(
            self,
            direct,
            method=method,
            update_along_time=update_along_time,
            tag=tag,
        )

    def load(self, direct, method="nc", tag=None):
        """
        Load function to load stored member variables.

        Parameters
        ----------
        direct : str
            directory at which the gcm_toolkit datasets are stored
        method : str, optional
            Should be the same method with which you stored the data
        tag: str, optional
            tag of the model that should be loaded.
            Will load all available models by default.
        """
        return raw.m_load(self, direct, method=method, tag=tag)

    # =============================================================
    #   Plotting Functions
    # =============================================================

    def isobaric_slice(self, var_key, pres, tag=None, **kwargs):
        """
        Plot an isobaric slice of the given quantity at the given pressure
        level. The user can specify the DataSet to be plotted by providing the
        corresponding tag.

        Parameters
        ----------
        var_key : str
            The key of the variable quantity that should be plotted.
        pres : float
            Pressure level for the isobaric slice to be plotted, expressed in
            the units specified in the dataset attributes (e.g., init of GCMT object).
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        # select the appropriate dataset
        dsi = self.get_one_model(tag)
        return gcmplt.isobaric_slice(dsi, var_key, pres, **kwargs)

    def time_evol(self, var_key, tag=None, **kwargs):
        """
        Function that plots the time evolution of a quantity
        in a 1D line collection plot, where the colorscale can be
        related to the time evolution.
        Note: var_key needs to contain data that is 2D in time and pressure.

        Parameters
        ----------
        dsi : DataSet
            A gcm_toolkit-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        dsi = self.get_one_model(tag)
        return gcmplt.time_evol(dsi, var_key, **kwargs)

    def zonal_mean(self, var_key, tag=None, **kwargs):
        """
        Plot a zonal mean average of a quantity for the given dataset.

        Parameters
        ----------
        dsi : DataSet
            A gcm_toolkit-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        tag : str, optional
            The tag of the dataset that should be plotted. If no tag is provided
            and multiple datasets are available, an error is raised.
        """
        # select the appropriate dataset
        dsi = self.get_one_model(tag)
        return gcmplt.zonal_mean(dsi, var_key, **kwargs)

    # ==============================================================================================
    #   Interfaces
    # ==============================================================================================

    def get_prt_interface(self, prt):
        """
        Constructs the interface to petitRADTRANS

        Parameters
        ----------
        prt: petitRADTRANS.Radtrans
            A fully initialized radtrans object to be used for the calculations

        Returns
        -------
        Interface: pRTInterface
            The interface object that is used to create phasecurves/spectra/etc
        """
        return PrtInterface(self, prt)

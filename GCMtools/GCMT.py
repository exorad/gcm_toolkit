# ==============================================================
#                       GCMtools Main Class
# ==============================================================
#  This class is the user interace of the GCMtools
#  functionalities. The goal is to have a clean and easy to use
#  environment for new users of GCMs while allowing direct
#  access to the data for more experienced users.
# ==============================================================

import xarray as xr
from collections import UserDict

from .read_in import m_read_from_mitgcm
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

    def __init__(self):
        """
        Constructor for the GCMtools class.

        Parameters
        ----------
        models : dict
            Dictionary containing all of the GCM datasets that have been read.
        """

        # Initialize empty dictionary to store all GCM models
        self.models = GCMDatasetCollection()

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
            All models in self.models, or only the one with the right tag.
        """
        # If no tag is given, return all models
        if tag is None:
            return self.models
        # If the tag is valid, return the corresponding model
        elif isinstance(tag, str):
            return self.models[tag]
        # If the tag is not a string, raise an error
        else:
            raise ValueError('The given tag is not a string.')


    def read_raw(self, gcm, data_path, iters='last', tag=None):
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
        tag : str
            Tag to reference the simulation in the collection of models.
        """

        # call the required GCM read-in method
        if gcm == 'MITgcm':
            ds = m_read_from_mitgcm(self, data_path, iters)
        else:
            raise ValueError('The selected GCM type "' + gcm +
                             '" is not supported')

        # if no tag is given, models are just numbered as they get added
        if tag is None:
            print('[WARN] -- No tag provided. This model is stored with tag: ' + str(len(self.models)))
            tag = str(len(self.models))

        # store tag in the dataset attributes
        ds.attrs['tag'] = tag

        # check if the dataset has all necessary GCMtools attributes
        if not is_the_data_basic(ds):
            raise ValueError('This dataset is not supported by GCMtools\n')

        # store dataset
        self.models[tag] = ds

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

        # if no tag is given, models are just numbered as they get added
        if tag is None:
            print('[WARN] -- No tag provided. This model is stored with tag: ' + str(len(self.models)))
            tag = str(len(self.models))

        # store tag in the dataset attributes    
        ds.attrs['tag'] = tag

        # check if the dataset has all necessary GCMtools attributes
        if not is_the_data_basic(ds):
            raise ValueError('This dataset is not supported by GCMtools\n')

        # store dataset
        self.models[tag] = ds

    def save(self):
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

        # Todo: add functionalities if required
        pass

    def load(self):
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

        # Todo: add functionalities if required
        pass

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

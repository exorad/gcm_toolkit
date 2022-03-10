# ==============================================================
#                       GCMtools Main Class
# ==============================================================
#  This class is the user interace of the GCMtools
#  functionalities. The goal is to have a clean and easy to use
#  environment for new users of GCMs while allowing direct
#  access to the data for more experienced users.
# ==============================================================

import xarray as xr

from .read_in import m_read_from_mitgcm
from .passport import is_the_data_basic

class GCMT:

    def __init__(self):
        """
        Constructor for the GCMtools class.

        Parameters
        ----------
        models : dict
            Dictionary containing all of the GCM datasets that have been read.

        Returns
        -------
        NoneType
            None
        """

        # Initialize empty dictionary to store all GCM models
        self.models = dict()
        pass

    def read_raw(self, gcm, data_path, key=None):
        """
        General read in function for GCM data

        Parameters
        ----------
        gcm : str
            Type of GCM, must be 'MITgcm'.
        data_path : str
            Folder path to the standard output of the GCM.
        key : str
            Key to reference the simulation in the collection of models.
        """

        if gcm == 'MITgcm':
            m_read_from_mitgcm(self, data_path)
        else:
            raise ValueError('The selected GCM type "' + gcm +
                             '" is not supported')

        # end read in function
        return

    def read_reduced(self, data_path, key=None):
        """
        Read in function for GCM data that has been reduced and saved according
        to the GCMtools format.

        Parameters
        ----------
        data_path : str
            Folder path to the reduced (GCMtools) data.
        key : str
            Key to reference the simulation in the collection of models.
        """
        # read dataset using xarray functionalities
        ds = xr.open_dataset(data_path)

        # check if the dataset has all necessary GCMtools attributes
        if not is_the_data_basic(ds):
            raise ValueError('The selected dataset is not supported by GCMtools\n' + gcm)

        # if no key is given, models are just numbered as they get added
        if key is None:
            key = len(models)

        # store dataset
        self.models[key] = ds

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

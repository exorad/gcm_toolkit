# ==============================================================
#                       GCMtools Main Class
# ==============================================================
#  This class is the user interace of the GCMtools
#  functionalities. The goal is to have a clean and easy to use
#  environment for new users of GCMs while allowing direct
#  access to the data for more experienced users.
# ==============================================================

from .read_in import m_read_from_mitgcm

class GCMT:

    def __init__(self):
        """
        Constructor for the GCMtools class.

        Parameters
        ----------
        example : int
            Short description of the variable.

        Returns
        -------
        NoneType
            None
        """

        # add stuff here
        pass

    def read_in(self, gcm, data_path):
        """
        General read in function for GCM data

        Parameters
        ----------
        gcm : str
            Type of GCM, must be 'MITgcm'.
        data_path : str
            Folder path to the standard output of the GCM.

        Returns
        -------
        NoneType
            None
        """

        if gcm == 'MITgcm':
            m_read_from_mitgcm(self, data_path)
        else:
            raise ValueError('The selected GCM type "' + gcm +
                             '" is not supported')

        # end read in function
        return

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

# ==============================================================
#                       GCMtools Plotting Class
# ==============================================================
#  This class incorporates all of GCMtools plotting
#  functionalities. We aim to have flexible and easy-to-use
#  plotting routines for the most common GCM data visualizations.
# ==============================================================

import matplotlib.pyplot as plt

class Plotting:
    """
    The GCMtools plotting class to collect all predefined
    plotting routines.

    Attributes
    ----------
    p_unit : str
        The pressure units (bar or pascal) used in plotting methods.

    Methods
    -------
    isobaric_slice:
        Plot the isobaric slice of a given dataset at a certain pressure level.
    zonal_mean:
        Plot the zonal mean quantity of a given dataset.
    """

    def __init__(self, p_unit='bar'):
        """
        Constructor for the GCMtools plotting class.

        Parameters
        ----------
        p_unit : str
            Sets the default pressure unit for all plotting routines to either
            'bar' (default) or 'Pa'.
        """
        # check if the given pressure unit is allowed (bar or pascal)
        if p_unit in ['bar', 'Pa']:
            self.p_unit = p_unit
        else:
            raise KeyError("[ERROR] The specified pressure unit (" + p_unit + ") is not allowed. Please use either 'bar' or 'Pa'.")

    def isobaric_slice(self, ds, var_key, p):
        """
        Plot an isobaric slice of the dataset and quantity at the given pressure
        level.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        p : float
            Pressure level for the isobaric slice to be plotted, expressed in
            the units specified in the init-method.
        """
        print('We are working in units: ' + self.p_unit)
        # TODO: implement
        pass

    def zona_mean(self, ds, var_key):
        """
        Plot a zonal mean average of a quantity for the given dataset.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
        var_key : str
            The key of the variable quantity that should be plotted.
        """
        # TODO: implement
        pass

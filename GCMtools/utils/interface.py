# ==============================================================
#                       GCMtools Interface Class
# ==============================================================
#  This class incorporates interfacing methods for GCMtools
#  datasets. The interfacing options are driven by the frequent
#  needs of users to transform GCM output to formats readable by
#  e.g.: - petitRADTRANS            (Molliere+2019)
#        - gCMCRT                   (Lee+2022)
#        - 1D and pseudo-2D PAC     (Agundez+2014)
#        - ...
# ==============================================================

class Interface:
    """
    The GCMtools interfacing class which implements common
    interfacing functionalities to different codes.

    Attributes
    ----------

    Methods
    -------
    """

    def to_2D_PAC(self, ds):
        """
        Transform the given GCM dataset to a pseudo-2D PAC
        input file.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
        """
        # TODO: implement later
        pass

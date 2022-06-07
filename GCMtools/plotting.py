# ==============================================================
#                       GCMtools Plotting Class
# ==============================================================
#  This class incorporates all of GCMtools plotting
#  functionalities. We aim to have flexible and easy-to-use
#  plotting routines for the most common GCM data visualizations.
# ==============================================================

import numpy as np
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

    def __init__(self):
        """
        Constructor for the GCMtools plotting class.

        Parameters
        ----------
        """


    def isobaric_slice(self, ds, var_key, p, time=-1, ax=None,
                       plot_windvectors=True, wind_kwargs=None, cbar_kwargs=None,
                       fs_labels=14, fs_ticks=11, title=None,
                       xlabel='Longitude (deg)', ylabel='Latitude (deg)',
                       **kwargs):
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
        if ax is None:
            fig= plt.figure()
            ax = plt.gca()
        if wind_kwargs is None:
            wind_kwargs = {}
        if cbar_kwargs is None:
            cbar_kwargs = {}

        # retrieve default units
        p_unit = ds.attrs.get('p_unit')
        time_unit = ds.attrs.get('time_unit')

        # By default, the plotted...
        #   - pressure is the nearest to p
        #   - iteration is the last one
        this_p = ds.Z.sel(Z=p, method="nearest").values
        this_time = ds.time.isel(time=-1).values
        ds2d = ds.isel(time=-1).sel(Z=p, method='nearest')

        # Simple plot (with xarray.plot.pcolormesh)
        plotted = ds2d[var_key].plot(add_colorbar=False, **kwargs)

        # make own colorbar, as the automatic colorbar is hard to customize
        cbar = plt.colorbar(plotted, ax=ax, **cbar_kwargs)
        if 'label' in cbar_kwargs:  # look if a colorbar label is given
            cbar_label = cbar_kwargs['label']
        else:
            cbar_label = var_key    # if not, use label of the plotted quantity
        cbar.set_label(cbar_label, fontsize=fs_labels)

        # Overplot the wind vectors if needed
        if plot_windvectors:
            self.plot_horizontal_wind(ds2d, ax=ax, **wind_kwargs)

        # set other plot qualities
        ax.set_aspect('equal')
        xt=np.arange(-180, 181, 45)
        yt=np.arange(-90, 91, 45)
        plt.xticks(ticks=xt, labels=[str(n)+r'$^\circ$' for n in xt], fontsize=fs_ticks)
        plt.yticks(ticks=yt, labels=[str(n)+r'$^\circ$' for n in yt], fontsize=fs_ticks)
        ax.set_xlabel(xlabel, fontsize=fs_labels)
        ax.set_ylabel(ylabel, fontsize=fs_labels)
        if title is None:
            title = f'p = {this_p:.2e} {p_unit}, time = {this_time:.0f} {time_unit}'
        ax.set_title(title, fontsize=fs_labels)


    def plot_horizontal_wind(self, ds, ax=None, sample_one_in=1, arrowColor='k'):
        """
        Plot the horizontal wind speeds as vector arrows.

        Parameters
        ----------
        ds : DataSet
            A GCMtools-compatible dataset where only latitude and longitude are
            non-singleton dimensions.

        Returns
        -------
        arrows : matplotlib.quiver.Quiver
            Quiver object that has been plotted.
        """
        if ax is None:
            fig= plt.figure()
            ax = plt.gca()

        # reduce the number of arrows plotted by sampling every one in n coordinates
        i = sample_one_in
        arrows = ax.quiver(ds.U.coords['lon'][::i], ds.U.coords['lat'][::i],
                           ds.U.values[::i, ::i], ds.V.values[::i, ::i],
                           pivot='mid', color=arrowColor)
        # TODO: perhaps use quiverkey to add a legend to the arrow length

        return arrows


    def zonal_mean(self, ds, var_key):
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

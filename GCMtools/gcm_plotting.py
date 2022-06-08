# ==============================================================
#                   GCMtools Plotting Library
# ==============================================================
#  This library incorporates all of GCMtools plotting
#  functionalities. We aim to have flexible and easy-to-use
#  plotting routines for the most common GCM data visualizations.
# ==============================================================

import numpy as np
import matplotlib.pyplot as plt

def isobaric_slice(ds, var_key, p, time=-1, lookup_method='exact', ax=None,
                   plot_windvectors=True, wind_kwargs=None, cbar_kwargs=None,
                   fs_labels=None, fs_ticks=None, title=None,
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
        the units specified in the dataset attributes (e.g., init of GCMT object).
    time : int, optional
        Timestamp that should be plotted. By default, the last time is
        selected.
    lookup_method : str, optional
        The look-up method that is used to slice along pressure:
        'exact' for exactly matching key look-up (default);
        'nearest' to pick out the nearest neighbour Z coordinate;
        'interpolate' for a linear interpolation along the Z axis.
    ax : matplotlib.axes.Axes, optional
        The axis on which you want your plot to appear.
    plot_windvectors : boolean, optional
        If True (default), plot the horizontal wind vectors over the colour
        mesh.
    wind_kwargs : dict, optional
        Additional keywords for the method plot_horizontal_wind.
    cbar_kwargs : dict, optional
        Additional keywords for the colorbar.
    fs_labels : int, optional
        Optionally set font size of the axis labels.
    fs_ticks : int, optional
        Optionally set font size of the tick labels.
    title : str, optional
        Title for the isobaric slice plot. By default, the selected pressure
        and time stamp of the slice are displayed.
    xlabel : str, optional
        X-axis label, longitude by default.
    ylabel : str, optional
        Y-axis label, latitude by default.
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

    # if no timestamp is given, pick the last available time
    if time == -1:
        time = ds.time.isel(time=-1).values
    # time-slice of the dataset
    # (note: the look-up method for time is always assumed to be exact)
    this_time = ds.time.sel(time=time).values
    ds = ds.sel(time=time)

    # isobaric slice based on the look-up method for pressure
    if lookup_method == 'exact':
        this_p = ds.Z.sel(Z=p).values
        ds2d = ds.sel(Z=p)
    elif lookup_method == 'nearest':
        this_p = ds.Z.sel(Z=p, method="nearest").values
        ds2d = ds.sel(Z=p, method='nearest')
    elif lookup_method == 'interpolate':
        this_p = ds.Z.interp(Z=p).values
        ds2d = ds.interp(Z=p)
    else:
        raise ValueError("Please enter 'exact', 'nearest', or 'interpolate' as Z lookup method.")

    # Simple plot (with xarray.plot.pcolormesh)
    plotted = ds2d[var_key].plot(add_colorbar=False, **kwargs)

    # make own colorbar, as the automatic colorbar is hard to customize
    cbar = plt.colorbar(plotted, ax=ax, **cbar_kwargs)
    cbar_label = cbar_kwargs.get('label', var_key)
    cbar.set_label(cbar_label, fontsize=fs_labels)

    # Overplot the wind vectors if needed
    if plot_windvectors:
        plot_horizontal_wind(ds2d, ax=ax, **wind_kwargs)

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


def plot_horizontal_wind(ds, ax=None, sample_one_in=1, arrowColor='k'):
    """
    Plot the horizontal wind speeds as vector arrows.

    Parameters
    ----------
    ds : DataSet
        A GCMtools-compatible dataset where only latitude and longitude are
        non-singleton dimensions.
    ax : matplotlib.axes.Axes, optional
        The axis on which you want your plot to appear.
    sample_one_in : int, optional
        If given, only one every nth coordinate will be sampled, to avoid
        overcrowding the figure.
    arrowColor : str, optional
        Specify the arrow color.

    Returns
    -------
    arrows : matplotlib.quiver.Quiver
        Quiver object that has been plotted.
    """
    if ax is None:
        fig= plt.figure()
        ax = plt.gca()

    # assert whether the dataset is 2d-horizontal, i.e. the vertical and time
    # dimensions only have a single coordinate
    if len(np.atleast_1d(ds.coords['Z'].values)) != 1:
        raise ValueError('This function only takes 2d-horizontal datasets. Select a single Z-coordinate first.')
    if len(np.atleast_1d(ds.coords['time'].values)) != 1:
        raise ValueError('This function only takes 2d-horizontal datasets. Select a single time-coordinate first.')

    # reduce the number of arrows plotted by sampling every one in n coordinates
    i = sample_one_in
    arrows = ax.quiver(ds.U.coords['lon'][::i], ds.U.coords['lat'][::i],
                       ds.U.values[::i, ::i], ds.V.values[::i, ::i],
                       pivot='mid', color=arrowColor)

    # TODO: perhaps use quiverkey to add a legend to the arrow length

    return arrows


def zonal_mean(ds, var_key):
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

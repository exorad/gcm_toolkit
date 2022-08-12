"""
==============================================================
                  gcm_toolkit Plotting Library
==============================================================
 This library incorporates all of gcm_toolkit plotting
 functionalities. We aim to have flexible and easy-to-use
 plotting routines for the most common GCM data visualizations.
==============================================================
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from ..core import writer as wrt
from ..core.const import VARNAMES as c


# pylint: disable=R0915,R0912
def isobaric_slice(
    dsi,
    var_key,
    pres,
    time=-1,
    lookup_method="exact",
    ax=None,
    plot_windvectors=True,
    wind_kwargs=None,
    cbar_kwargs=None,
    add_colorbar=True,
    fs_labels=None,
    fs_ticks=None,
    title=None,
    xlabel="Longitude (deg)",
    ylabel="Latitude (deg)",
    contourf=False,
    **kwargs,
):
    """
    Plot an isobaric slice of the dataset and quantity at the given pressure
    level.

    Parameters
    ----------
    dsi : DataSet
        A gcm_toolkit-compatible dataset of a 3D climate simulation.
    var_key : str
        The key of the variable quantity that should be plotted.
    pres : float
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
    add_colorbar: bool, optional,
        Optionally decide if you want a colorbar or don't
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
    contourf: bool, optional
        Decide if you want to do a contourplot or a pcolormesh plot
    """

    # print information
    wrt.write_status("STAT", "Plot Isobaric slices")
    wrt.write_status("INFO", "Variable to be plotted: " + var_key)
    wrt.write_status("INFO", "Pressure level: " + str(pres))

    if ax is None:
        plt.figure()
        ax = plt.gca()
    if wind_kwargs is None:
        wind_kwargs = {}
    if cbar_kwargs is None:
        cbar_kwargs = {}

    font_labels = dict(fontsize=fs_labels) if fs_labels is not None else {}
    font_ticks = dict(fontsize=fs_ticks) if fs_ticks is not None else {}

    # retrieve default units
    p_unit = dsi.attrs.get("p_unit")
    time_unit = dsi.attrs.get("time_unit")

    # if no timestamp is given, pick the last available time
    if time == -1:
        time = dsi[c["time"]].isel(**{c["time"]: -1}).values
    # time-slice of the dataset
    # (note: the look-up method for time is always assumed to be exact)
    dsi = dsi.sel(**{c["time"]: time})

    # isobaric slice based on the look-up method for pressure
    if lookup_method == "exact":
        this_p = dsi[c["Z"]].sel(**{c["Z"]: pres}).values
        ds2d = dsi.sel(**{c["Z"]: pres})
    elif lookup_method == "nearest":
        this_p = dsi[c["Z"]].sel(**{c["Z"]: pres}, method="nearest").values
        ds2d = dsi.sel(**{c["Z"]: pres}, method="nearest")
    elif lookup_method == "interpolate":
        this_p = dsi[c["Z"]].interp(**{c["Z"]: pres}).values
        ds2d = dsi.interp(**{c["Z"]: pres})
    else:
        raise ValueError(
            "Please enter 'exact', 'nearest', or 'interpolate' as Z lookup"
            " method."
        )

    # Simple plot (with xarray.plot.pcolormesh)
    if contourf:
        plotted = ds2d[var_key].plot.contourf(
            add_colorbar=False, ax=ax, **kwargs
        )
    else:
        plotted = ds2d[var_key].plot.pcolormesh(
            add_colorbar=False, ax=ax, **kwargs
        )

    # make own colorbar, as the automatic colorbar is hard to customize
    if add_colorbar:
        cbar = plt.colorbar(plotted, ax=ax, **cbar_kwargs)
        cbar_label = cbar_kwargs.get("label", var_key)
        cbar.set_label(cbar_label, **font_labels)

    # Overplot the wind vectors if needed
    if plot_windvectors:
        plot_horizontal_wind(ds2d, ax=ax, **wind_kwargs)

    # set other plot qualities
    if not hasattr(ax, "projection"):
        ax.set_aspect("equal")
        xtp = np.arange(-180, 181, 45)
        ytp = np.arange(-90, 91, 45)
        ax.set_xticks(xtp)
        ax.set_yticks(ytp)
        ax.set_xticklabels([str(n) + r"$^\circ$" for n in xtp], **font_ticks)
        ax.set_yticklabels([str(n) + r"$^\circ$" for n in ytp], **font_ticks)

    ax.set_xlabel(xlabel, **font_labels)
    ax.set_ylabel(ylabel, **font_labels)
    if title is None:
        if time_unit == "iter":
            # need to convert time from nanosecond like datatype to iters
            time_string = f"{1e-9 * float(time):.0f}"
        else:
            time_string = f"{time}"
        title = f"p = {this_p:.2e} {p_unit}, time = {time_string} {time_unit}"
    ax.set_title(title, **font_labels)


def plot_horizontal_wind(
    dsi, ax=None, sample_one_in=1, arrow_color="k", windstream=False, **kwargs
):
    """
    Plot the horizontal wind speeds as vector arrows.

    Parameters
    ----------
    dsi : DataSet
        A gcm_toolkit-compatible dataset where only latitude and longitude are
        non-singleton dimensions.
    ax : matplotlib.axes.Axes, optional
        The axis on which you want your plot to appear.
    sample_one_in : int, optional
        If given, only one every nth coordinate will be sampled, to avoid
        overcrowding the figure.
    arrowColor : str, optional
        Specify the arrow color.
    windstream: bool, optional
        Specify if you want to plot a streamfunction or arrows.

    Returns
    -------
    arrows : matplotlib.quiver.Quiver
        Quiver object that has been plotted.
    """

    # print information
    wrt.write_status("STAT", "Plot horizontal winds")

    if ax is None:
        plt.figure()
        ax = plt.gca()

    # assert whether the dataset is 2d-horizontal, i.e. the vertical and time
    # dimensions only have a single coordinate
    if len(np.atleast_1d(dsi.coords[c["Z"]].values)) != 1:
        raise ValueError(
            "This function only takes 2d-horizontal datasets. Select a single"
            " Z-coordinate first."
        )
    if len(np.atleast_1d(dsi.coords[c["time"]].values)) != 1:
        raise ValueError(
            "This function only takes 2d-horizontal datasets. Select a single"
            " time-coordinate first."
        )

    if not windstream or not hasattr(ax, "projection"):
        # reduce the number of arrows plotted by sampling every one in n coordinates
        i = sample_one_in
        arrows = ax.quiver(
            dsi[c["U"]].coords[c["lon"]][::i],
            dsi[c["U"]].coords[c["lat"]][::i],
            dsi[c["U"]].values[::i, ::i],
            dsi[c["V"]].values[::i, ::i],
            pivot="mid",
            color=arrow_color,
            **kwargs,
        )

        # perhaps use quiverkey to add a legend to the arrow length
    else:
        u_data, v_data = dsi[c["U"]].values, dsi[c["V"]].values
        speed = np.sqrt(u_data**2 + v_data**2)
        lwt = (speed / speed.max()) ** 0.5

        arrows = ax.streamplot(
            dsi[c["U"]].coords[c["lon"]],
            dsi[c["U"]].coords[c["lat"]],
            u_data,
            v_data,
            linewidth=lwt,
            color=arrow_color,
            density=sample_one_in,
            **kwargs,
        )

    return arrows


def _multiline(xsi, ysi, c_map, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xsi : iterable container of x coordinates
    ysi : iterable container of y coordinates
    c_map : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xsi, ysi)]
    lcn = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lcn.set_array(np.asarray(c_map))

    # add lines to axes and rescale
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lcn)
    ax.autoscale()
    return lcn


def time_evol(
    dsi,
    var_key,
    ax=None,
    fs_labels=None,
    cbar_kwargs=None,
    add_colorbar=True,
    title=None,
    xlabel=None,
    ylabel="Z",
    add_ylabel_unit=True,
    **kwargs,
):
    """
    Function that plots the time evolution of a quantity
    in a 1D line collection plot, where the colorscale can
    be related to the time evolution.
    Note: var_key needs to contain data that is 2D in time and pressure.

    Parameters
    ----------
    dsi : DataSet
        A gcm_toolkit-compatible dataset of a 3D climate simulation.
    var_key : str
        The key of the variable quantity that should be plotted.
    ax : matplotlib.axes.Axes, optional
        The axis on which you want your plot to appear.
    fs_labels : int, optional
        Optionally set font size of the axis labels.
    cbar_kwargs : dict, optional
        Additional keywords for the colorbar.
    add_colorbar: bool, optional
        Decide if you want a colorbar
    title : str, optional
        Title for the isobaric slice plot. By default, the selected pressure
        and time stamp of the slice are displayed.
    xlabel: str, optional
        Label for x
    ylabel: str, optional
        Label for y
    add_ylabel_unit: bool, optional
        Optionally decide, if you want to add a unit to ylabel.

    Returns
    -------
    l: LineCollection
        the collection of plotted lines

    """

    # print information
    wrt.write_status("STAT", "Plot horizontal winds")
    wrt.write_status("INFO", "Variable to be plotted: " + var_key)

    p_unit = dsi.attrs.get("p_unit")
    time_unit = dsi.attrs.get("time_unit")

    if cbar_kwargs is None:
        cbar_kwargs = {}

    font_labels = dict(fontsize=fs_labels) if fs_labels is not None else {}

    if ax is None:
        plt.figure()
        ax = plt.gca()

    # Test the dimension of the quantity.
    # Need dsi to be one two-dimensional (time + pressure)
    if len(dsi[var_key].dims) != 2:
        raise ValueError(
            "Time evolution plots only work with two dimensional data (time +"
            " pressure)!"
        )

    xsv, ysv = [], []
    for t_it in dsi.time:
        x_it = dsi[var_key].sel(**{c["time"]: t_it})
        y_it = dsi[x_it.dims[0]]
        xsv.append(x_it)
        ysv.append(y_it)

    l_out = _multiline(xsi=xsv, ysi=ysv, c_map=dsi.time, ax=ax, **kwargs)

    # make own colorbar, as the automatic colorbar is hard to customize
    if add_colorbar:
        cbar = plt.colorbar(l_out, ax=ax, **cbar_kwargs)
        cbar_label = cbar_kwargs.get("label", f"time ({time_unit})")
        cbar.set_label(cbar_label, **font_labels)

    # set other plot qualities
    if title is None:
        title = f"timeevolution of {var_key}"
    ax.set_title(title, **font_labels)

    if add_ylabel_unit:
        ylabel = ylabel + f" ({p_unit})"
    if xlabel is None:
        xlabel = var_key

    ax.set_xlabel(xlabel, **font_labels)
    ax.set_ylabel(ylabel, **font_labels)

    # Invert y-axis and set scale to log
    ax.set_yscale("log")
    ax.invert_yaxis()
    return l_out


def zonal_mean(
    dsi,
    var_key,
    time=-1,
    ax=None,
    cbar_kwargs=None,
    fs_labels=None,
    xlabel="Latitude (deg)",
    ylabel="Z",
    add_ylabel_unit=True,
    title=None,
    add_colorbar=True,
    contourf=False,
    **kwargs,
):
    """
    Plot a zonal mean average of a quantity for the given dataset.

    Parameters
    ----------
    dsi : DataSet
        A gcm_toolkit-compatible dataset of a 3D climate simulation.
    var_key : str
        The key of the variable quantity that should be plotted.
    time : int, optional
        Timestamp that should be plotted. By default, the last time is
        selected.
    ax : matplotlib.axes.Axes, optional
        The axis on which you want your plot to appear.
    cbar_kwargs : dict, optional
        Additional keywords for the colorbar.
    fs_labels : int, optional
        Optionally set font size of the axis labels.
    xlabel: str, optional
        Label for x
    ylabel: str, optional
        Label for y
    add_ylabel_unit: bool, optional
        Optionally decide, if you want to add a unit to ylabel.
    title : str, optional
        Title for the isobaric slice plot. By default, the selected pressure
        and time stamp of the slice are displayed.
    add_colorbar: bool, optional
        Optionally decide if you want a colorbar or don't
    contourf: bool, optional
        Decide if you want to do a contourplot or a pcolormesh plot

    """

    # print information
    wrt.write_status("STAT", "Plot zonal mean")
    wrt.write_status("INFO", "Variable to be plotted: " + var_key)

    if ax is None:
        plt.figure()
        ax = plt.gca()
    if cbar_kwargs is None:
        cbar_kwargs = {}

    font_labels = dict(fontsize=fs_labels) if fs_labels is not None else {}

    # retrieve default units
    p_unit = dsi.attrs.get("p_unit")
    time_unit = dsi.attrs.get("time_unit")

    # if no timestamp is given, pick the last available time
    if time == -1:
        time = dsi[c["time"]].isel(**{c["time"]: -1}).values
    # time-slice of the dataset
    # (note: the look-up method for time is always assumed to be exact)
    # this_time = time
    zmean = dsi[var_key].sel(**{c["time"]: time}).mean(dim=c["lon"])

    # Simple plot (with xarray.plot.pcolormesh)
    if contourf:
        plotted = zmean.plot.contourf(
            add_colorbar=False, ax=ax, x=c["lat"], **kwargs
        )
    else:
        plotted = zmean.plot.pcolormesh(
            add_colorbar=False, ax=ax, x=c["lat"], **kwargs
        )

    # make own colorbar, as the automatic colorbar is hard to customize
    if add_colorbar:
        cbar = plt.colorbar(plotted, ax=ax, **cbar_kwargs)
        cbar_label = cbar_kwargs.get("label", var_key)
        cbar.set_label(cbar_label, **font_labels)

    # set other plot qualities
    if title is None:
        if time_unit == "iter":
            # need to convert time from nanosecond like datatype to iters
            time_string = f"{1e-9 * float(time):.0f}"
        else:
            time_string = f"{time}"

        title = f"time = {time_string} {time_unit}"
    ax.set_title(title, **font_labels)

    if add_ylabel_unit:
        ylabel = ylabel + f" ({p_unit})"

    ax.set_xlabel(xlabel, **font_labels)
    ax.set_ylabel(ylabel, **font_labels)

    # Invert y-axis and set scale to log
    ax.set_yscale("log")
    ax.invert_yaxis()

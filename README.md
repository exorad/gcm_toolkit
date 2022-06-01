# GCMtools
> *GCMtools* is an open-source python package to read, post-process, and plot 3D GCM data. The goal is to have an easy to use invironment for new GCM users, while allowing for direct access to the data for more experienced users.

*GCMtools* originated from the need for a consistent and easy-to-share methodology for reading and analyzing the 3D circulation models of exoplanet atmospheres. Currently, the focus is on data produced by the [MITgcm](http://mitgcm.org/) and its implementations for exoplanets (see [Carone et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.3582C/abstract) and [Schneider et al., 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220209183S/abstract)). But we envision this package as a general tool for the analysis of data from different 3D GCMs.

*GCMtools* is currently under development. Feel free to contribute!

## Features (Current/Planned)

* Multiple read-in methods (currently only *MITgcm*)
* 'Tag' system for easy referencing of multiple models in memory
* Making use of [xarray](https://docs.xarray.dev/en/stable/)'s intuitive system of multi-D arrays
* Quick plotting routines
    * isobaric slices
    * zonal means
    * evolution over time
* Direct data access for custom plotting and post-processing

## Main Authors
* [Aaron Schneider](https://github.com/AaronDavidSchneider)
* [Sven Kiefer](https://github.com/Kiefersv)
* [Robin Baeyens](https://github.com/robinbaeyens)

## Get Started

### Installation

`git clone https://github.com/exorad/GCMtools`  

`pip install -e <directory containing setup.py>`

### Packages needed

xarray, netcdf4, xmitgcm

cubedsphere (see https://cubedsphere.readthedocs.io/en/latest/index.html)

### Running the first example

Even if you do not have GCM output data matching one of the supported models, you can still get started with *GCMtools* by loading pre-reduced datasets (in NetCDF format). Have a look at the example dataset `example_dataset.nc` and plot it with `python get_started.py`.

## Miscellaneous Documentation

* All quantities are expressed in _SI units_.
* Latitudes run from -90&deg; to 90&deg; ; Longitudes run from -180&deg; to 180&deg;.
* The substellar point for tidally locked planets is located at (0&deg;, 0&deg;).
* The radius of the planet is defined at the bottom boundary of the vertical grid.

## Contact

Feel free to contact one of the main authors: **Aaron Schneider**, **Sven Kiefer**, or **[Robin Baeyens](https://www.linkedin.com/in/robin-baeyens/)**.
Project url: https://github.com/exorad/GCMtools

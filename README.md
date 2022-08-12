# gcm_toolkit

[![codecov](https://codecov.io/gh/exorad/gcm_toolkit/branch/main/graph/badge.svg?token=OTEPCHABXI)](https://codecov.io/gh/exorad/gcm_toolkit)
[![PyPi version](https://badgen.net/pypi/v/pip/)](https://pypi.org/project/gcm-toolkit/)
[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](http://gcm-toolkit.readthedocs.io/?badge=latest)

> *gcm_toolkit* is an open-source python package to read, post-process, and plot 3D GCM data. The goal is to have an
> easy to use invironment for new GCM users, while allowing for direct access to the data for more experienced users.

*gcm_toolkit* originated from the need for a consistent and easy-to-share methodology for reading and analyzing the 3D
circulation models of exoplanet atmospheres. Currently, the focus is on data produced by
the [MITgcm](http://mitgcm.org/) and its implementations for exoplanets (
see [Carone et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.3582C/abstract)
and [Schneider et al., 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220209183S/abstract)). But we envision this
package as a general tool for the analysis of data from different 3D GCMs.

*gcm_toolkit* is currently under development. Feel free to contribute!

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

Read the [docs](https://gcm-toolkit.readthedocs.io/en/latest/)

### Installation

`pip install gcm-toolkit`

(Note the dash)

### Packages needed

Please install the following packages beforehand:

**general:**
xarray, netcdf4

**exPERT/MITgcm:**
xmitgcm, xgcm, f90nml,
[cubedsphere](https://cubedsphere.readthedocs.io/en/latest/index.html)

### Running the first example

Check the get started in the [docs](https://gcm-toolkit.readthedocs.io/en/latest/notebooks/demo.html)

## Miscellaneous Documentation

* All quantities are expressed in _SI units_.
* Latitudes run from -90&deg; to 90&deg; ; Longitudes run from -180&deg; to 180&deg;.
* The substellar point for tidally locked planets is located at (0&deg;, 0&deg;).
* The radius of the planet is defined at the bottom boundary of the vertical grid.

## Contact

Feel free to contact one of the main authors: **[Aaron Schneider](https://aaron-schneider.de)**
, **[Sven Kiefer](https://kiefersv.github.io/)**, or **[Robin Baeyens](https://www.linkedin.com/in/robin-baeyens/)**.

Project url: [https://github.com/exorad/gcm_toolkit](https://github.com/exorad/gcm_toolkit)

.. GCMTools documentation master file, created by
   sphinx-quickstart on Mon Jan 11 14:18:49 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GCMTools's documentation!
=======================================
Library for post processing of GCM data

Philosophy
----------
``GCMTools`` aims at simplifying the post-processing of GCM data.
It is mainly built with the purpose of generalizing the process of dealing with exoplanet GCM data.

We appreciate that every GCM comes with different output formats and structures.
However, we find that most users are interested in similar diagnostics, such as plotting isobaric slices or plotting zonal means.

``GCMTools`` is therefore tailored to abstract away the differences in the output of GCM data.

Capabilities
------------

.. warning::
   ``GCMTools`` is under active development. If you would like to see your GCM supported by ``GCMTools``, please contact one of the authors.

**Implemented functionality:**

- Plotting (isobaric slices, zonal means, etc.)
- readin capabilities for exorad
- load and save capabilities of the ``GCMTools`` abstraction
- interface to ``petitRADTRANS`` for phasecurves, etc.

**Planned functionality:**

- interface to ``gCMCRT`` for complex radiative transfer calculations

.. toctree::
   :maxdepth: 1
   :caption: Contents

   installation
   usage
   cli
   notebooks/demo.ipynb
   notebooks/petitRADTRANS.ipynb

   
Contributers
------------
- `Robin Baeyens <https://github.com/robinbaeyens/>`_
- `Aaron Schneider <https://github.com/aarondavidschneider/>`_
- `Sven Kiefer <https://github.com/kiefersv/>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Usage
=====

The functionalities of ``gcm_toolkit`` is likely to expand in the future, since ``gcm_toolkit`` is in active development.
However, ``gcm_toolkit`` can already be used to perform a variety of tasks.
These docs explain how ``gcm_toolkit`` can be used to load data and do plots with the data.

GCMT class
----------
The GCMT class is the main point of user interaction.
It is used to load data into ``xarray`` datasets, which can then be used for analysis.

.. code-block:: python

        from gcm_toolkit import GCMT

        # Load data
        tools = GCMT()
        tools.read_raw(..., tag='model_name')
        # tools.read_raw(..., tag='model_name_2')  # if you have multiple models
        ds = tools['model_name']

The GCMT class can also be used for quick plotting of data, since it wraps a few plotting functions, which are outlined in :ref:`Plotting`.

.. note::
   
   You can also load netcdf data:
   
   .. code-block:: python
   
        tools.read_reduced(..., tag='model_name')


.. autoclass:: gcm_toolkit.GCMT
    :members: __init__, get, get_models, models, read_raw, read_reduced, load, save


Postprocessing
--------------
The GCMT class has quite a few methods that can calculate diagnostics such as total energy, location of the RCB and so on.
You find their documentation here.

Example for the total energy:

.. code-block:: python

        from gcm_toolkit import GCMT
        import gcm_toolkit.gcm_plotting as gcmp

        # Load data
        tools = GCMT()
        tools.read_raw(..., tag='model_name')
        E = tools.add_total_energy(var_key_out='E', tag='model_name')

        # You will also have the energy stored in the dataset if var_key_out is not None:
        ds = tools['model_name']
        E == ds.E

.. autoclass:: gcm_toolkit.GCMT
    :members: add_theta, add_horizontal_average, add_total_energy, add_meridional_overturning, add_rcb, add_total_momentum


Plotting
--------
``gcm_toolkit`` includes some very simple (but yet powerful and highly customizable) plotting functions that can be used in two different ways.
You can either use the ``gcm_toolkit`` class to invoke these functions or call these functions directly.

Example:

.. code-block:: python

        from gcm_toolkit import GCMT
        import gcm_toolkit.gcm_plotting as gcmp

        # Load data
        tools = GCMT()
        tools.read_raw(..., tag='model_name')
        ds = tools['model_name']

        # Plot data with GCMT
        tools.isobaric_slice(tag = 'model_name', p=1e-2, lookup_method='nearest', var_key='T', wind_kwargs={'windstream':False, 'sample_one_in':2})

        # Plot data with GCMTtools.gcm_plotting
        tools.isobaric_slice(ds, p=1e-2, lookup_method='nearest', var_key='T', wind_kwargs={'windstream':False, 'sample_one_in':2})

Both ways are equally powerful and flexible.

.. autofunction:: gcm_toolkit.gcm_plotting.isobaric_slice

.. autofunction:: gcm_toolkit.gcm_plotting.zonal_mean

.. autofunction:: gcm_toolkit.gcm_plotting.time_evol


Interface
---------
``gcm_toolkit`` comes with some interfaces to other codes.
We currently have support for ``petitRADTRANS`` to calculate spectra and phasecurves with ``prt_phasecure``.

.. autoclass:: gcm_toolkit.utils.interface.PrtInterface
    :members: __init__, set_data, chem_from_poorman, calc_phase_spectrum, phase_curve



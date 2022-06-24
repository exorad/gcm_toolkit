Usage
=====

The functionalities of ``GCMTools`` is likely to expand in the future, since ``GCMTools`` is in active development.
However, ``GCMTools`` can already be used to perform a variety of tasks.
These docs explain how ``GCMTools`` can be used to load data and do plots with the data.

GCMT class
----------
The GCMT class is the main point of user interaction.
It is used to load data into ``xarray`` datasets, which can then be used for analysis.

.. code-block:: python

        from GCMTtools import GCMT

        # Load data
        gcmt = GCMT()
        gcmt.read_raw(..., tag='model_name')
        # gcmt.read_raw(..., tag='model_name_2')  # if you have multiple models
        ds = gcmt.get_models('model_name')

The GCMT class can also be used for quick plotting of data, since it wraps a few plotting functions, which are outlined in :ref:`Plotting`.

.. autoclass:: GCMtools.GCMT
    :members: __init__, get_models, models, add_horizontal_average, read_raw, read_reduced, load, save


Plotting
--------
``GCMTools`` includes some very simple (but yet powerful and highly customizable) plotting functions that can be used in two different ways.
You can either use the ``GCMT`` class to invoke these functions or call these functions directly.

Example:

.. code-block:: python

        from GCMTtools import GCMT
        import GCMTtools.gcm_plotting as gcmp

        # Load data
        gcmt = GCMT()
        gcmt.read_raw(..., tag='model_name')
        ds = gcmt.get_models('model_name')

        # Plot data with GCMT
        gcmt.isobaric_slice(tag = 'model_name', p=1e-2, lookup_method='nearest', var_key='T', wind_kwargs={'windstream':False, 'sample_one_in':2})

        # Plot data with GCMTtools.gcm_plotting
        gcmp.isobaric_slice(ds, p=1e-2, lookup_method='nearest', var_key='T', wind_kwargs={'windstream':False, 'sample_one_in':2})

Both ways are equally powerful and flexible.

.. autofunction:: GCMtools.gcm_plotting.isobaric_slice

.. autofunction:: GCMtools.gcm_plotting.zonal_mean

.. autofunction:: GCMtools.gcm_plotting.time_evol

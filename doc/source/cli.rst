Command Line Interface
======================

``gcm_toolkit`` also provides some functionality to be used from the command line.
The functionalites explained below are scripts, which are automatically installed into your environment, when installing ``gcm_toolkit``.
You can use them in your terminal by hitting the name of the script (e.g.:).
Information on the required arguments can be found below, or by using the ``-h`` flag.


convert_to_gcmt
---------------

The ``convert_to_gcmt`` script is meant to enable an intermediary step between postprocessing/dataanalysis/plotting and simulation outputs.
``convert_to_gcmt`` reads the raw input of the GCM output and saves it to an intermediary dataformat (netcdf or zarr).

The script could be used in the following way:

1. Do climate simulation (remotely on cluster)
2. Run ``convert_to_gcmt`` to convert the data to netcdf
3. Download the netcdf data and load it into ``GCMT``
4. Work on the data

**Usage:**

.. code-block:: bash

    $ convert_to_gcmt -h
    usage: convert_to_gcmt [-h] [-c CONFIG]

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            specify the path to the config file that holds the
                            config for converting to xarray

The use of a config file (yaml) for the script is optional, but highly recommended, since script will fallback to defaults, which might or might not work for your setup.
If you do not provide any name for a config file, the script will search for ``convert.yaml`` in your current working directory.

The following options are available in the ``convert.yaml`` file:

.. code-block:: yaml

    iters: "all"                 # list of iterations to load
    prefix: ["T","U","V","W"]    # The file prefixes that you want to read
    tag: "nc_convert"            # GCMT tag that should be created for your conversion
    gcm: "MITgcm"                # The GCM for which we want to load the raw data
    load_existing: False         # If we want to check for already converted input and just extend if available
    data_path: "run"             # Path to the raw data
    save_path: "results"         # Path at which the converted data should be stored
    method: "nc"                 # Output format used for the conversion
    update_along_time: False     # Only relevant if method == "zarr". Checkout gcm_toolkit.GCMT.save for more info
    # (anything else to be passed to gcm_toolkit.GCMT.read_raw)

.. Note:: All of the other arguments are input for :meth:`gcm_toolkit.GCMT.read_raw`

Checkout :meth:`gcm_toolkit.GCMT.read_raw`, :meth:`gcm_toolkit.GCMT.load` and :meth:`gcm_toolkit.GCMT.save` to understand the usage of the above parameters.
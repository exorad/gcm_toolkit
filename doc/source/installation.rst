Installation
============
There are several ways to install the package.
We recommend using anaconda.
``gcm_toolkit`` has been tested to work on linux and mac.

.. note::
    Using Windows is not recommended and will not work in combination with MITgcm raw data.

Using anaconda
----------------
Create an anaconda environment.

.. code-block:: bash

    conda create -n 'GCM' -y
    conda activate GCM

Install a few dependencies (optional include the following packages: ``jupyterlab`` and ``cartopy``.)

.. code-block:: bash

    conda install -c conda-forge cubedsphere jupyterlab cartopy -y
    pip install git+https://github.com/MITgcm/xmitgcm.git


Install ``gcm_toolkit``.

.. code-block:: bash

    pip install gcm-toolkit

.. note::
    If you prefer to install a development version, you can use the following command:

    .. code-block:: bash

        git clone https://github.com/exorad/gcm_toolkit.git
        pip install -e gcm_toolkit

    We recommend to install ``pre-commit`` into your git hooks.
    This will automatically format the code to meet common style guidelines.
    Its as easy as:

    .. code-block:: bash

        cd gcm_toolkit
        pre-commit install

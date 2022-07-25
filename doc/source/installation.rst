Installation
============
There are several ways to install the package.
We recommend using anaconda.
``gcmt`` has been tested to work on linux and mac.

.. note::
    Using Windows is not recommended and will not work in combination with MITgcm raw data.

Using anaconda
----------------
Create an anaconda environment.

.. code-block:: bash

    conda create -n 'MITGCM' -y
    conda activate MITGCM

Install a few dependencies (optional include the following packages: ``jupyterlab`` and ``cartopy``.)

.. code-block:: bash

    conda install -c conda-forge cubedsphere jupyterlab cartopy -y
    pip install git+https://github.com/MITgcm/xmitgcm.git


Install ``gcmt``.

.. code-block:: bash

    pip install gcmt

.. note::
    If you prefer to install a development version, you can use the following command:

    .. code-block:: bash

        git clone https://github.com/exorad/gcmt.git
        pip install -e gcmt

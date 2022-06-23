Installation
============

Using anaconda:

.. code-block:: bash

    conda create -n 'MITGCM' -y
    conda activate MITGCM
    conda install -c conda-forge cubedsphere jupyterlab cartopy -y
    pip install git+https://github.com/MITgcm/xmitgcm.git
    pip install git+https://github.com/exorad/GCMtools.git

This installation already includes the following packages: ``jupyterlab`` and ``cartopy``.


.. _oneapi-installation:

*************************
Intel oneAPI Installation
*************************

This installation guide works on Fedora Linux.

Add Intel oneAPI Repository
===========================

First, import Intel's GPG key and add the official repository:

.. code-block:: bash

    # Import Intel GPG key
    sudo rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

    # Create Intel oneAPI repository configuration file
    echo -e "[intel-oneapi]\nname=Intel(R) oneAPI Repository\nbaseurl=https://yum.repos.intel.com/oneapi\nenabled=1\ngpgcheck=1\ngpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" | sudo tee /etc/yum.repos.d/intel-oneapi.repo

Install Intel oneAPI Components
===============================

Intel MKL Math Kernel Library:

.. code-block:: bash

    sudo dnf install intel-oneapi-mkl intel-oneapi-mkl-devel

Intel MPI Library:

.. code-block:: bash

    sudo dnf install intel-oneapi-mpi intel-oneapi-mpi-devel

Set Environment Variables
=========================

Edit ``~/.bashrc`` or ``~/.zshrc`` file and add:

.. code-block:: bash

    source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

Uninstallation
==============

.. code-block:: bash

    # Uninstall specific component
    sudo dnf remove intel-oneapi-mkl

    # Uninstall all oneAPI components
    sudo dnf remove intel-oneapi-\*

Notes
=====

- Intel oneAPI requires a relatively new kernel version, Fedora 36 or later is recommended
- Installation may take several minutes to over ten minutes depending on network speed
- Full installation may require approximately 10-20GB of disk space
- It is recommended to update the system with ``sudo dnf update`` before installation
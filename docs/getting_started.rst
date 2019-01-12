.. _getting_started:


***************
Getting started
***************

.. _installing-docdir:

Third-party package dependence
===============================

  * **MPICH** --- an MPI implementation library, available at http://www-unix.mcs.anl.gov/mpi/mpich

  * **GSL** --- the GNU Scientific Library, downloaded at http://www.gnu.org/software/gs

  * **LAPACKE** --- the C-interface of LAPACK, downloaded at http://www.netlib.org/lapack/

  * **CDNest** --- Diffusive nested sampling, downloaded at https://github.com/LiyrAstroph/CDNest

Note that in Linux system, there are package managers that can install the above libraries convienently. If so, use them. In this case, the libraries usually are installed in standard environment path. Otherwise, any of the above libraries is not installed in standard locations on your system, the ``Makefile`` provided with the code may need slight adjustments.

Compiling
=============================

Edit the configurations in ``Makefile`` to be consistent with your system if necessary. Then compile the package with the command

.. code:: bash

   make

This creates an executable file ``mica2``.

Running
=============================

To run the package in a parallel computer/cluster, use the following command: 

.. code:: bash

   mpiexec -n np ./mica2 param

where ``np`` is the number of cores and ``param`` is the paramter file.

Parameter file
=============================

A typical parameter file looks like::

  #
  # lines starting with "#" are regarded as comments and are neglected
  #

  FileDir                   ./
  DataFile                  data/data.txt
 
  FlagUniformVarParams      0                # whether each dataset has the same variability parameters

  FlagUniformTranFuns       0                # whether each dataset has the same line parameters.
                                             # note that different lines have different parameters.

  LagLimitLow               0.0              # lower limit of the range of time lag to be explored
  LagLimitUpp               10.0             # upper limit of the range of time lag to be explored

  NumGaussianLow            2                # lower limit of number of Gaussians
  NumGaussianUpp            2                # upper limit of number of Gaussians

  FlagConSysErr             1                # 0, not include systematic error of continuum; 1, include
  FlagLineSysErr            1                # 0, not include systematic error of line; 1, include

Data format
==============================

``mica2`` reads data files with a format as::

  # 1
  # 171:269
  56690.6100           3.4270          0.0640     % continuum, 171 lines
  56691.5400           3.5450          0.0650
  ...
  56864.8600           4.3310          0.0740
  56865.9200           4.7080          0.0780

  56698.3570           2.1900          0.0560     % line, 269 lines
  56699.5590           2.2000          0.0580
  ...
  56830.1490           2.3000          0.0650
  56830.4200           2.2900          0.0660

The first line starting with "#" specifies the number of datasets. Here one dataset contains one continuum light curve (the driving source) and several line light curves (at least one). The second line starting with "#" specifies the numbers of points in light curves of continuum and lines, which are separated by ":". 

The next follows data of light curves, going by datasets. For each dataset, the first block is continuum light curve and then line light curves successively. In each data block, the three columns are time, flux, and error, respectively. Blocks/datasets are separated by a blank line.

In the above example, there is one dataset and it contains 171 points in continuum light curve and 269 point in one line light curve. If your data have 2 datasets,  the first dataset has 2 line light curves while the second dataset has one light curve, the data file should be formated as::
  
  # 2
  # 171:130:90
  # 150:122
  56690.6100           3.4270          0.0640     % continuum of 1st dataset, 171 lines
  56691.5400           3.5450          0.0650
  ...
  56864.8600           4.3310          0.0740
  56865.9200           4.7080          0.0780

  56698.3570           2.1900          0.0560     % 1st line of 1st dataset, 130 lines
  56699.5590           2.2000          0.0580
  ...
  56830.1490           2.3000          0.0650
  56830.4200           2.2900          0.0660

  56698.3570           2.1900          0.0560     % 2nd line of 1st dataset, 90 lines
  56699.5590           2.2000          0.0580
  ...
  56830.1490           2.3000          0.0650
  56830.4200           2.2900          0.0660

  56690.6100           3.4270          0.0640     % continuum of 2nd dataset, 150 lines
  56691.5400           3.5450          0.0650
  ...
  56864.8600           4.3310          0.0740
  56865.9200           4.7080          0.0780

  56698.3570           2.1900          0.0560     % line of 2nd dataset, 122 lines
  56699.5590           2.2000          0.0580
  ...
  56830.1490           2.3000          0.0650
  56830.4200           2.2900          0.0660

As you can see, the numbers of lines in each datasets do not needs to be the same.
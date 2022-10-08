.. _getting_started:


***************
Getting Started
***************

.. _installing-docdir:

Third-party package dependence
===============================

  * **MPICH** --- an MPI implementation library, available at http://www-unix.mcs.anl.gov/mpi/mpich

  * **GSL** --- the GNU Scientific Library, downloaded at http://www.gnu.org/software/gs

  * **LAPACKE** --- the C-interface of LAPACK, downloaded at http://www.netlib.org/lapack/

  * **CDNest** --- Diffusive nested sampling, downloaded at https://github.com/LiyrAstroph/CDNest

Note that in Linux system, there are package managers that can install the above libraries convienently (except CDNest). 
If so, use them. In this case, the libraries usually are installed in standard environment path. Otherwise, any of the above 
libraries is not installed in standard locations on your system, the ``Makefile`` provided with the code may need slight adjustments.

Compiling
=============================

Edit the configurations in ``Makefile`` to be consistent with your system if necessary. Then compile the package with the command

.. code:: bash

   make

This creates an executable file ``mica2``.

Running
=============================

First create two subdirectories ``data/`` and ``param`` in the current working directory. All the output files will be placed 
into ``data/``. The subdirectory ``param`` is used to place options for CDNest. 

To run the package in a parallel computer/cluster, use the following command: 

.. code:: bash

   mpiexec -n np ./mica2 param/param

where ``np`` is the number of cores and ``param`` is the paramter file, stored in the directory ``param/``.
This will also generate CDNest option files ``OPTIONSCON`` and ``OPTIONS1D`` in the subdirectory ``param/``.

If the results are not as good as expected, one may want to modify options for Markov-chain Monte Carlo sampling.
There are two ways. The first way is directly editing the parameter file (such as ``param/param`` in the above; see below).
The second way is editing the above generated option file ``OPTIONS1D`` and transfer it to ``mica2`` in the command line as

.. code:: bash

  mpiexec -n np ./mica2 param/param param/OPTIONS1D

where ``OPTIONS1D`` is an options file stored in the directory ``param/``, 
see :ref:`cdnest_label` for the detail.

Parameter file
=============================

A typical parameter file looks like::

  #
  # lines starting with "#" are regarded as comments and are neglected
  # if want to turn on the line, remove the beginning "#"
  
  #==============================================================
  
  FileDir                   ./
  DataFile                  data/sim_data.txt
  
  TypeTF                    0                # 0: Gaussian
                                             # 1: Top-hat
  
  MaxNumberSaves            2000             # number of MCMC sampling steps
  
  FlagUniformVarParams      0                # whether each dataset has the same variability parameters
  
  FlagUniformTranFuns       0                # whether each dataset has the same line parameters.
                                             # note that different lines have different parameters.
  
  FlagLongtermTrend         0                # Longterm trend in light curves, use a polynomial to fit 
                                             # input the order of the polynomial, e.g.,
                                             # 0, constant  (default)
                                             # 1, linear line 
                                             # 2, conic line
                                             # Use the default if you do not know this.
  
  LagLimitLow               0.0              # lower limit of the range of time lag to be explored
  LagLimitUpp               100.0              # upper limit of the range of time lag to be explored
                                             # can be negative
  
  NumCompLow                2                # lower limit of number of Gaussians/tophats
  NumCompUpp                2                # upper limit of number of Gaussians/tophats
  
  FlagConSysErr             0                # 0, not include systematic error of continuum; 1, include
  FlagLineSysErr            0                # 0, not include systematic error of line; 1, include
  
  TypeLagPrior              0                # type of lag prior for each Gaussians/tophats.
                                             # 0,  limit0 < lag0 < lag1 < lag2 <... < limit1
                                             #
                                             # 1,  limit0 + 0*width < lag0 < limit0 + 1*width
                                             #     limit0 + 1*width < lag1 < limit0 + 2*width
                                             #     ...
                                             #     width = (limit1 - limit0)/num_comp
                                             #
                                             # 2,  lags fixed at specific values, no limit on Guassian sigma/tophat width
                                             #     lag0 = limit0 + 0*dlag
                                             #     lag1 = limit0 + 1*dlag
                                             #     ...
                                             #     dlag = (limit1 - limit0)/(num_comp-1)
                                             #     
                                             # 3,  lags fixed at specific values
                                             #     Gaussian sigma ranges at (dlag/2, dlag), tophat wdith=dlag/2
                                             #     lag0 = limit0 + 0*dlag
                                             #     lag1 = limit0 + 1*dlag
                                             #     ...
                                             #     dlag = (limit1 - limit0)/(num_comp-1)
                                             #     better to set a large  mumber of components
  
  
  #==================================================================
  # options for cdnest sampling
  # use the default values or do not turn thme on IF NOT familiar with them
  
  # PTol                     0.1             # likelihood tolerance in loge
  # NumberParticles          1               # number of particles
  # NewLevelIntervalFactor   2               # new level interval
  # SaveIntervalFactor       2               # particular saving interval
  # ThreadStepsFactor        10              # thread steps before communications between cores
  # MaxNumberLevels          0               # maximum number of levels; unlimited for 0
  # BacktrackingLength       10.0            # backforward tracking length (lambda)
  # StrengthEqualPush        100.0           # strength to force equal push (beta)
  
  #===================================================================

see :ref:`cdnest_label` for the detail of CDNest options.

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

The next follows data of light curves, going by datasets. For each dataset, the first block is continuum light curve and then line light curves successively. In each data block, the three columns are time, flux, and error, respectively. **Blocks/datasets are separated by a blank line.**

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

Output
=================================
``mica2`` outputs the following main files in the folder ``data/``:  
  
  * **posterior_sample1d.txt_xx**

    posterior sample for parameters. The postfix "_xx" means the number of Gaussians. 
    The order of parameters in posterior sample file is arranged as: 

    (systematic error of continuum, sigmad, taud) * number of datasets

    (systematic error of line, (gaussian amplitude, center, sigma) * number of gaussians * number of line datasets) * number of datasets

    sigmad, taud, gaussian amplitude and sigma are in logarithm scale; systematic errors (x) are dimensionless, defined as  x = log(1+err/err_data), where err is the real systematic error and err_data is the mean measurement error of the data.

  * **pall.txt_xx**

    reconstruction of datasets, with the same format as the input data.
  
  * **pline.txt_xx_compyy** (applicable with ``-d`` option)
    
    decomposed light curves for each Gaussian component, with the same format as the input data. **yy** (a number) indicates which Gaussian component.
    Note that the continuum light curve is not decomposed and only line light curves are decomposed. 

  * **para_names_line.txt_xx**
    
    parameters and their priors.

  * **evidence.txt**
    
    Bayesian evidence for each number of Gaussians explored.

In the end of running, ``mica2`` prints the obtained Bayesian evidence for each number of Gausssians explored.

Plotting
========================
There is a Python script **plotfig.py** provided in the package that can be used to plot the results. 
Run it with 

.. code:: python

  ptyhon plotfig.py param/param

This will generate a PDF file **fig_xx.pdf** in the subdirectory **data/**.
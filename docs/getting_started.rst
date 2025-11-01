.. _getting_started:


***************
Getting Started
***************

.. _installing-docdir:

Third-party package dependence
===============================

  * **MPICH** --- an MPI implementation library, available at http://www-unix.mcs.anl.gov/mpi/mpich. 
    
    Note that in some cases, the package **hwloc** is not automatically installed when installing mpich. One needs to install it 
    and its development package **hwloc-devel** to use mpi.

    There is another popular MPI library `OpenMPI <https://www.open-mpi.org/>`_. If one sticks to OpenMPI, make sure to use the same 
    MPI library between compiling and running MICA.

  * **GSL** --- the GNU Scientific Library, downloaded at http://www.gnu.org/software/gs

  * **LAPACKE** --- the C-interface of LAPACK, downloaded at http://www.netlib.org/lapack/

  * **CMake** ---- the compilation tool, downloaded at https://cmake.org/download/

Note that one does not need to compile the above packages from the source code. Using package managers will be more convenient.

- On Linux system, there are package managers that can install the above libraries convienently. 
  If so, use them, e.g., for Fedora/Redhat distribution 

  .. code-block:: bash 

    dnf install mpich mpich-devel hwloc hwloc-devel gsl gsl-devel lapack lapack-dvel

  In this case, the libraries usually are installed in standard environment path. Otherwise, any of the above 
  libraries is not installed in standard locations in your system, the compiling configurations below may need slight adjustments.

- On MacOS, one can use `homebrew <https://brew.sh>`_ package manager to install the required packages

  .. code-block:: bash 

    brew install mpich hwloc gsl lapack

  To make lapack and hwloc findable by pkgconfig, export the paths of ``lapack.pc`` and ``hwloc.pc`` to the 
  environment  ``PKG_CONFIG_PATH``. For example, if lapack is installed at ``/opt/homebrew/opt/lapack`` and 
  hwloc is installed at ``/opt/homebrew/opt/hwloc``, then execute the commands

  .. code-block:: bash

    export PKG_CONFIG_PATH=/opt/homebrew/opt/lapack/lib/pkgconfig:$PKG_CONFIG_PATH
    export PKG_CONFIG_PATH=/opt/homebrew/opt/hwloc/lib/pkgconfig:$PKG_CONFIG_PATH
  
  or put the above commands into ``.bashrc``.  If the above commands do not work, one can also add the paths to 
  the environments ``LDFLAGS`` and ``CFLAGS`` in the terminal 

  .. code-block:: bash 
    
    export LDFLAGS=-L/opt/homebrew/opt/lapack/lib:$LDFLAGS
    export CPPFLAGS=-I/opt/homebrew/opt/lapack/include:$CPPFLAGS
    
    export LDFLAGS=-L/opt/homebrew/opt/hwloc/lib:$LDFLAGS
    export CPPFLAGS=-I/opt/homebrew/opt/hwloc/include:$CPPFLAGS

For Intel CPUs, the Intel OneAPI MKL library provides optimized interfaces to LAPACKE and BLAS libraries. Using MKL libraray can improve
the running speed and it is therefore highly recommended. MICA automatically checks the system environment variable ``MKLROOT`` to 
determine whether MKL libraray has been installed. If yes, MICA will by default use the interfaces of the MKL library; if not, the original 
LAPACKE and BLAS libraraies will be used. See `Intel OneAPI MKL <https://www.intel.com/content/www/us/en/docs/onemkl/get-started-guide/2025-0/overview.html>`_
for the installation. Please keep in mind that one needs to source the MKL's variable in `bashrc`, so as to let it callable from MICA.

Compiling
=============================

Exectuable Binary Version: ``mica2``
------------------------------------

Edit CMake compiling configurations to be consistent with your system if necessary, using the command

.. code:: bash

  ccmake .

This will draw out a CMake GUI and edit the configurations accordingly. **If this is your first installation, the GUI might be 
empty and type ``c`` to do initial configurating.** Please refer to 
`CMake <https://cmake.org/cmake/help/book/mastering-cmake/index.html>`_ for more details about the useage of CMake GUI.

After this step, compile the package with the command

.. code:: bash

   make

This creates an executable file ``mica2``.

If your system does not have lastest CMake, you can use Makefile provided in the package to do compiling. 
First edit the configurations in ``Makefile_old`` to be consistent with your system's setting, and then execute 
the command 

.. code:: bash

   make -f Makefile_old


Python Callable Version: ``pymica``
-----------------------------------

The Python package `mpi4py <https://pypi.org/project/mpi4py/>`_ is required. Install it using ``pip``

.. code-block:: bash 

  python -m pip install mpi4py

Note that ``pip`` keeps previously built wheel files on its cache for future reuse. 
If you want to reinstall the mpi4py package using a different or updated MPI implementation, 
you have to either first remove the cached wheel file with

.. code-block:: bash 
  
  python -m pip cache remove mpi4py

or ask ``pip`` to disable the cache:

.. code-block:: bash
  
  python -m pip install --no-cache-dir mpi4py

Now install MICA using the command 

.. code:: bash 

  python setup.py install --user 

This will generate a Python package ``pymica`` and install it to the user's Python package sites. 
In the folder ``tests/python``, the Python script ``example.py`` shows how to use pymica.

Running with Binary Version
=============================

First create two subdirectories ``data/`` and ``param`` in the current working directory. All the output files will be placed 
into ``data/``. The subdirectory ``param`` is used to place options for CDNest. 

To run the package in a parallel computer/cluster, use the following command, e.g.: 

.. code:: bash

   mpiexec -n 6 ./mica2 param/param  # here use 6 cores, change it to the numbers you want

where ``param`` is the paramter file, stored in the directory ``param/``.
This will also generate CDNest option files ``OPTIONSCON`` and ``OPTIONS1D`` in the subdirectory ``param/``.

If the results are not as good as expected, one may want to modify options for Markov-chain Monte Carlo sampling.
There are two ways. The first way is directly editing the parameter file (such as ``param/param`` in the above; see below).
The second way is editing the above generated option file ``OPTIONS1D`` and transfer it to ``mica2`` in the command line as

.. code:: bash

  mpiexec -n 6 ./mica2 param/param param/OPTIONS1D # here use 6 cores, change it to the numbers you want

where ``OPTIONS1D`` is an options file stored in the directory ``param/``, 
see :ref:`cdnest_label` for the detail.

Parameter File
=============================

A typical parameter file looks like::

  #
  # lines starting with "#" are regarded as comments and are neglected
  # if want to turn on the line, remove the beginning "#"
  # note that some options are optinal
  
  #==============================================================
  
  FileDir                   ./
  DataFile                  data/IRAS_year5.txt

  TypeModel                 0                # 0: general model
                                             # 1: pmap, photometric RM
                                             # 2: vmap, use a virtual driving light curve
                                             # 3: mmap, mixture of TF types
                                             # default: 0

  TypeTF                    0                # 0: Gaussian
                                             # 1: Top-hat
                                             # 2: Gamma function (k=2)
                                             # 3: Exponential
                                             # default: 0
  
  StrTypeTFMix              30               # string for the TF types of mixture mode (mmap)
                                             # e.g., "02": gaussian-gamma; "20": gamma-gaussian
                                             # valid if TypeModel == 3 (mmap)

  MaxNumberSaves            1000             # number of MCMC sampling steps
                                             # default: 2000

  FlagUniformVarParams      0                # whether each dataset has the same variability parameters
                                             # default: 0

  FlagUniformTranFuns       0                # whether each dataset has the same line parameters.
                                             # note that different lines have different parameters.
                                             # default: 0

  FlagLongtermTrend         0                # Longterm trend in light curves, use a polynomial to fit 
                                             # input the order of the polynomial, e.g.,
                                             # 0, constant  (default)
                                             # 1, linear line 
                                             # 2, conic line
                                             # Use the default if you do not know this.

  LagLimitLow               0.0              # lower limit of the range of time lag to be explored
  LagLimitUpp               300.0            # upper limit of the range of time lag to be explored
                                             # can be negative

  #WidthLimitLow            1.0              # lower and upper limit of lag width
  #WidthLimitUpp            50.0             # by default, MICA determines the limits automatically.
                                             # if unsatifactory, turn on these options.

  FlagLagPositivity         0                # whether force Gaussians overall located at non-negative lags
                                             # 0: no;  1: yes
                                             # default: 0

  FlagNegativeResp          1                # whether turn on negative response
                                             # 0, no; 1, yes
                                             # default: 0

  NumCompLow                2                # lower limit of number of Gaussians/tophats
  NumCompUpp                2                # upper limit of number of Gaussians/tophats
                                             # default: 1, 1

  FlagConSysErr             0                # 0, not include systematic error of continuum; 1, include
  FlagLineSysErr            1                # 0, not include systematic error of line; 1, include
                                             # defaul: 0, 0
  
  NumPointRec              200               # number of points in reconstruction for each light curve
                                             # note: a too large number causes very slow reconstruction
                                             # default: 200
  
  TimeRecLowExt            0                 # extend time range of reconstruction
  TimeRecUppExt            0                 # if detnote the original time range as [t1, t2]
                                             # then the options lead to [t1+TimeRecLowExt, t2+TimeRecUppExt]
                                             # default: 0, 0
  
  #StrWidthPrior           [1:10:5:20]       # width priors if the default priors not good enough
                                             # format: [width1_1:width1_2:width2_1:width2_2...]
                                             # "WidthLimitLow" and "WidthLimitUpp" no longer applicable 
                                             # default: None

  TypeLagPrior              1                # type of lag prior for each Gaussians/tophats.
                                             # default: 0
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
                                             #
                                             # 4,  user specified with "StrLagPrior"

  StrLagPrior             [0:10:10:50]       # valid if TypeLagPrior==4
                                             # format: [lag1_1:lag1_2:lag2_1:lag2_2...]
                                             # "LagLimitLow" and "LagLimitUpp" no longer applicable

  StrRatioPrior           [1.0e-3:1.0]       # the response ratio of 2nd to 1st component
                                             # valid if TypeModel == 1 (pmap)
                                             # format: [ratio_1:ratio_2]
                                             # default: [1.0e-3:1.0]
  
  FlagGap                 0                  # whether include seasonal gap
                                             # 0: no; 1: yes.
                                             # default: 0

  #StrGapPrior            [182.6:140.0]      # gap priors if the default priors are not good enough
                                             # valid when FlagGap == 1
                                             # format: [gap_center_set1:gap_width_set1:gap_center_set2:gap_width_set2...]
                                             # gap_center_set1: gap center for 1st dataset (+n*year will also be included)
                                             # gap_width_set1:  gap width for 1st dataset
                                             # default: None
  #==================================================================
  # options for cdnest sampling
  # use the default values or do not turn them on IF NOT familiar with them
  
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

Running with Python Version
==============================
In Python environment, import mica and other necessary packages as, 

.. code:: python

  from mpi4py import MPI
  import numpy as np
  import matplotlib.pyplot as plt
  import pymica

Then initialize MPI environment as 

.. code:: python

  # initiate MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

If one did not create the formated data file (see below), one could directly load the light curves 
and feed them to MICA as  

.. code:: python

  if rank == 0:
    con = np.loadtxt("cont.txt")
    line= np.loadtxt("line.txt")

    # make a data dict 
    data_input = {"set1":[con, line]}

    # if multiple datasets, e.g., 
    #data_input = {"set1":[con1, line1], "set2":[con2, line2]}

    # if a dataset has multiple lines, e.g.,
    #data_input = {"set1":[con, line1, line2]}
  else:
    data_input = None 

  data_input = comm.bcast(data_input, root=0)

  #create a model
  #there are two ways
  #1) one way from the param file

  #model = pymica.gmodel(param_file="param/param_input")

  #2) the ohter way is through the setup function
  
  # type: gmodel(), pmap(), vmap()
  model = pymica.gmodel()

  # type： gaussian, tophat, gamma, exp
  model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[1, 1], max_num_saves=2000)

If one already has created the formatted data file (see blow), one can directly input the file name as 

.. code:: python

  model.setup(data_file="file_name", type_tf='gaussian', lag_limit=[0, 100], number_component=[1, 1], max_num_saves=2000)


After the above initialization, run the code as 

.. code:: python
  
  #run mica
  model.run()

  #posterior run, only re-generate posterior samples, do not run MCMC
  # model.post_run()

  #do decomposition for the cases of multiple components 
  # model.decompose()

  # plot results
  if rank == 0:
    
    # plot results, doshow controls whether showing the results on screen
    # 
    model.plot_results(doshow=True, tf_lag_range=None, hist_lag_range=None, show_pmax=True) 
    model.post_process()  # generate plots for the properties of MCMC sampling 

    # get the full sample 
    # sample is a list, each element contains an array of posterior samples
    # sample[0] is for the case of number_component[0]
    # sample[1] is for the case of number_component[1] 
    # ...
    sample = model.get_posterior_sample()

See :ref:`python_label` for a detailed guideline.

Data Format
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

.. code:: bash

  ptyhon plotfig.py --param param/param

This will generate a PDF file **fig_xx.pdf** in the subdirectory **data/**.  Use the following command 
to print help information about this script.  

.. code:: bash

  python plotfig.py --help

Testing
========================
To test ``mica2``, change to the subdirectory ``tests/`` and there are several suites of tests to guide the useage of ``mica2``. 
See :ref:`tests_label` for more details.
MICA2
=====
.. image:: https://readthedocs.org/projects/mica2/badge/?version=latest
  :target: https://mica2.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Read `the documentation <http://mica2.readthedocs.io/>`_.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you make use of this code, please cite our paper `Li, Y.-R., Wang, J.-M., & Bai, J.-M. 2016, ApJ, 831, 206 <http://adsabs.harvard.edu/abs/2016ApJ...831..206L>`_.

+++++++++++++++++
A Quick Tutorial
+++++++++++++++++

Types of Transfer Functions
---------------------------
MICA2 supports serveral types of tranfer functions:

- Gaussian 
- Tophat 
- Gamma distribution (k=2)

Exectuable Binary Version: mica2
--------------------------------

Prepare the parameter file (e.g., named **param**) and input data file (placed in sub-folder **data/**), 
and then run mica in a terminal , e.g., 

.. code-block:: bash 

  mpiexec -n 6 ./mica2 param 

Here **-n 6** means using 6 cores. Change the number according to your need.

Python Module Version: pymica
--------------------------------

Edit a Python script named, e.g., example.py, as the following.

.. code-block:: Python

  from mpi4py import MPI
  import numpy as np
  import pymica
  import matplotlib.pyplot as plt
  
  # initiate MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  
  # load data
  if rank == 0:
    data = np.loadtxt("./sim_data.txt")
    con = data[:126, :]
    line= data[126:, :]
    
    # make a data dict 
    data_input = {"set1":[con, line]}
  
    # if multiple datasets, e.g., 
    #data_input = {"set1":[con1, line1], "set2":[con2, line2]}
    
    # if a dataset has multiple lines, e.g.,
    #data_input = {"set1":[con, line1, line2]}
  else:
    data_input = None 
  
  data_input = comm.bcast(data_input, root=0)
  
  model = pymica.gmodel()
  model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[2, 2], max_num_saves=200)
  # if using tophats, set type_tf='tophat'
  
  #the full arguments are 
  #model.setup(data_file=None, data=None,
  #            type_tf='gaussian', max_num_saves=2000, 
  #            flag_uniform_var_params=False, flag_uniform_tranfuns=False,
  #            flag_trend=0, flag_lag_posivity=False,
  #            lag_limit=[0, 100], number_component=[1, 1],
  #            flag_con_sys_err=False, flag_line_sys_err=False,
  #            type_lag_prior=0)
  
  #run mica
  model.run()
  
  # plot results
  if rank == 0:
    # get the full sample 
    # sample is a list, each element contains an array of posterior samples
    # sample[0] is for the case of number_component[0]
    # sample[1] is for the case of number_component[1] 
    # ...
    sample = model.get_posterior_sample()
  
    # get the posterior sample of time lags of the "line" in the dataset "set"
    # timelag is a list, each element contains an array of posterior samples
    # timelag[0] is for the case of number_component[0]
    # timelag[1] is for the case of number_component[1]
    # ...
    timelag = model.get_posterior_sample_timelag(set=0, line=0) 
    plt.plot(timelag[0][:, 0])
    plt.plot(timelag[0][:, 1])
    plt.show()

    # get the posterior sample of widths of the "line" in the dataset "set"
    # width is a list, each element contains an array of posterior samples
    # width[0] is for the case of number_component[0]
    # width[1] is for the case of number_component[1]
    # ...
    width = model.get_posterior_sample_width(set=0, line=0)  
    plt.plot(width[0][:, 0])
    plt.plot(width[0][:, 1])
    plt.show() 
  
    model.plot_results() # plot results
    model.post_process()  # generate plots for the properties of MCMC sampling 
  

Run this script using the terminal command as 

.. code-block:: bash

  mpiexec -n 6 python example.py 

If you want to use only one core, just run as 

.. code-block:: bash

  python example.py 

Photometric Reverberation Mapping
---------------------------------

MICA2 can also do reverberation mapping analysis between two photometric light curves, in which 
the photometric bands may contain broad-line emissions or other components so that there may 
exist multiple responses. For simplicity, MICA2 assumes that the driving photometric light curve 
does not contain those contaminations and purely reflects continuum variations.

.. code-block:: python

  from mpi4py import MPI
  import numpy as np
  import pymica
  import matplotlib.pyplot as plt
  
  # initiate MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  
  # load data
  if rank == 0:
    data = np.loadtxt("./sim_data.txt")
    band1 = data[:126, :]
    band2 = data[126:, :]
    
    # make a data dict 
    data_input = {"set1":[band1, band2]}
  else:
    data_input = None 
  
  data_input = comm.bcast(data_input, root=0)
    
  model = pymica.pmap()
  model.setup(data=data_input, type_tf='gaussian', max_num_saves=2000, lag_prior=[[-5, 5],[0, 50]], ratio_prior=[0.01, 0.5])
  # if using tophats, set type_tf='tophat'
    
  #run mica
  model.run()
  
  # plot results
  if rank == 0:
     
    model.plot_results() # plot results
    model.post_process()  # generate plots for the properties of MCMC sampling


Virtual Reverberation Mapping
-----------------------------

MICA2 also provides a ``vmap`` mode to do reverberation mapping analysis with a virtual driving light curve. This mode applies 
in cases where the dirving light curve cannot be chosen or the driving light curve has a poor qaulity that is not suitable to act 
as the dirving one.

To this end, MICA2 assumes that the virtual drving light curve follows a DRW process with a variation amplitude (:math:`\sigma`) of 0.1 
and has a time lag of zero with respect to the first light curve of the input data. The remaining analysis is trival and 
similar to the normal modes.

.. code-block:: python
  
  from mpi4py import MPI
  import numpy as np
  import pymica
  import matplotlib.pyplot as plt
  
  # initiate MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  # load data
  if rank == 0:
    lc0 = np.empty(0)  # virtual light curve, empty
    lc1 = np.loadtxt("g.txt")
    lc2 = np.loadtxt("r.txt")
    
    # make a data dict 
    data_input = {"set1":[lc0, lc1, lc2]}
  else:
    data_input = None 
  
  data_input = comm.bcast(data_input, root=0)
  
  model = pymica.vmap()
  model.setup(data=data_input, type_tf='gaussian', lag_limit=[-2, 5], number_component=[1, 1], max_num_saves=2000)
  # if using tophats, set type_tf='tophat'
  # see the documentation for the format of vmap data.

  #run mica
  model.run()
  
  # plot results
  if rank == 0:
    
    model.plot_results() # plot results
    model.post_process()  # generate plots for the properties of MCMC sampling 
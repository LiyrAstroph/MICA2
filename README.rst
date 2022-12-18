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

Exectuable Binary Version: mica2
--------------------------------

Prepare the parameter file (e.g., named **param**) and input data file (placed in sub-folder **data/**), 
and then run mica in a terminal , e.g., 

.. code-block:: bash 

  mpiexec -n 6 ./mica2 param 

Here **-n 6** means using 6 cores.

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
  data = np.loadtxt("./sim_data.txt")
  con = data[:126, :]
  line= data[126:, :]
  
  # make a data dict 
  data_input = {"set1":[con, line]}
  
  model = pymica.model()
  model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[2, 2], max_num_saves=200)
  
  #the full arguments are 
  #model.setup(data_file=None, data=None,
  #            type_tf='gaussian', max_num_saves=2000, 
  #            flag_uniform_var_params=False, flag_uniform_tranfuns=False,
  #            flag_trend=0, flag_lag_posivity=False,
  #            lag_limit=[0, 100], number_component=[1, 1],
  #            flag_con_sys_err=False, flag_line_sys_err=False,
  #            type_lag_prior=0): 
  
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
  
    # get the posterior sample of time lags of the line in dataset 
    # timelag is a list, each element contains an array of posterior samples
    # timelag[0] is for the case of number_component[0]
    # timelag[1] is for the case of number_component[1]
    # ...
    timelag = model.get_posterior_sample_timelag(set=0, line=0)  
  
    model.plot_results() # plot results
    model.postprocess()  # generate plots for the properties of MCMC sampling 
  

Run this script using the terminal command as 

.. code-block:: bash

  mpiexec -n 6 python example.py 

If you want to use only one core, just run as 

.. code-block:: bash

  python example.py 
.. _python_label:

Running with Python Version
===========================

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

If one already has created the formatted data file (see :ref:`getting_started`), one can directly input the file name as 

.. code:: python

  model.setup(data_file="file_name", type_tf='gaussian', lag_limit=[0, 100], number_component=[1, 1], max_num_saves=2000)

The full arguments for setup function looks like 

.. code:: python

  model.setup(data_file=None, data=None,
              type_tf='gaussian', max_num_saves=2000, 
              flag_uniform_var_params=False, flag_uniform_tranfuns=False,
              flag_trend=0, flag_lag_posivity=False,
              lag_limit=[0, 100], number_component=[1, 1],
              width_limit=[0.1, 100],
              nd_rec=200, trec_ext=[0, 0],
              flag_con_sys_err=False, flag_line_sys_err=False,
              type_lag_prior=0, lag_prior=[[0, 50]],
              width_prior[[1, 40]],
              num_particles=1, thread_steps_factor=1, 
              new_level_interval_factor=1, save_interval_factor=1,
              lam=10, beta=100, ptol=0.1, 
              max_num_levels=0)

The meaning of these arguements are the same as in the binary version (see :ref:`getting_started`). 
Most of the arguments are optional and are not must be specified.

After the above initialization, execuate the following steps 

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

Up to this point, the run is completed and a serious of output files are generated in `data/` 
directory (see :ref:`getting_started`).  A PDF file "fig_xx.pdf" is created to plot the reconstructed 
light curves and the corresponding transfer functions. A PDF file "cdnest_xx.pdf" is created to show
the statistics of MCMC sampling.

One then directly read in the posterior sample of parameters from the file "posterior_sample1d.txt_xx".
Alternatively, one can also load the posterior sample using 

.. code:: python
  
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

  # get the posterior sample of widths of the "line" in the dataset "set"
  # width is a list, each element contains an array of posterior samples
  # width[0] is for the case of number_component[0]
  # width[1] is for the case of number_component[1]
  # ...
  width = model.get_posterior_sample_width(set=0, line=0) 

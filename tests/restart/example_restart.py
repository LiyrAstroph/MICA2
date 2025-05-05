# 
# an exemplary script to show how to resume from the last run
#
# include two scripts,
# 1) example_initial.py: do a run 
# 2) example_restart.py: resume from the run

# this is the script for resuming the last run
# you need to firstly execute "example_initial.py" to do a run
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
# 

from mpi4py import MPI
import numpy as np
import pymica
import matplotlib.pyplot as plt
import shutil

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

#create a model
#there are two ways
#1) one way from the param file

#model = pymica.gmodel(param_file="param/param_input")

#2) the ohter way is through the setup function

model = pymica.gmodel()
# use Gaussians
model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[1, 2], max_num_saves=1000)
# gaussian:    type_tf = "gaussian"
# tophat:      type_tf = "tophat"
# exponential: type_tf = "exp"
# gamma (k=2): type_tf = "gamma"

#the full arguments are 
#model.setup(data_file=None, data=None,
#            type_tf='gaussian', max_num_saves=2000, 
#            flag_uniform_var_params=False, flag_uniform_tranfuns=False,
#            flag_trend=0, flag_lag_posivity=False,
#            lag_limit=[0, 100], number_component=[1, 1],
#            width_limit=[0.1, 100],
#            nd_rec=200,
#            flag_con_sys_err=False, flag_line_sys_err=False,
#            type_lag_prior=0, lag_prior=[[0, 50]],
#            width_prior[[1, 40]],
#            num_particles=1, thread_steps_factor=1, 
#            new_level_interval_factor=1, save_interval_factor=1,
#            lam=10, beta=100, ptol=0.1, 
#            max_num_levels=0)


# prepare the snapshorts
# this will copy "restart_dnest1d.txt_xx_yyyy" to "restart_dnest1d.txt_xx",
#                xx: number of components; yyyy: steps from which to resume
#
model.prepare(last_steps=500)

#resume from last run
model.restart()

# plot results
if rank == 0:
  
  # plot results, doshow controls whether showing the results on screen
  # 
  model.plot_results(doshow=True, tf_lag_range=None, hist_lag_range=None, hist_bins=None) 
  model.post_process()  # generate plots for the properties of MCMC sampling 

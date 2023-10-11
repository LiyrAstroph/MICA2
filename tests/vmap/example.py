# 
# an examplary script to show the useage of pymica
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
# 

from mpi4py import MPI
import numpy as np
import pymica
import matplotlib.pyplot as plt

# initiate MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# load data
if rank == 0:
  lc0 = np.empty(0)
  lc1 = np.loadtxt("mcg08_g.txt")
  lc2 = np.loadtxt("mcg08_r.txt")
  
  # make a data dict 
  data_input = {"set1":[lc0, lc1, lc2]}
else:
  data_input = None 

data_input = comm.bcast(data_input, root=0)

#create a model
#there are two ways
#1) one way from the param file

#model = pymica.gmodel(param_file="param/param_input")

#2) the ohter way is through the setup function

model = pymica.vmap()

# use Gaussians
model.setup(data=data_input, type_tf='gaussian', lag_limit=[-2, 5], number_component=[1, 1], max_num_saves=5000)

# or use tophats
#model.setup(data=data_input, type_tf='tophat', lag_limit=[0, 100], number_component=[2, 2], max_num_saves=2000)

#the full arguments are 
#model.setup(data_file=None, data=None,
#            type_tf='gaussian', max_num_saves=2000, 
#            flag_uniform_var_params=False, flag_uniform_tranfuns=False,
#            flag_trend=0, flag_lag_posivity=False,
#            lag_limit=[0, 100], number_component=[1, 1],
#            width_limit=[0.1, 100],
#            flag_con_sys_err=False, flag_line_sys_err=False,
#            type_lag_prior=0, lag_prior=[[0, 50]],
#            num_particles=2, thread_steps_factor=2, 
#            new_level_interval_factor=2, save_interval_factor=2,
#            lam=10, beta=100, ptol=0.1, 
#            max_num_levels=0)

#run mica
model.run()

#posterior run, only re-generate posterior samples, do not run MCMC
# model.post_run()

#do decomposition for the cases of multiple components 
#model.decompose()

# plot results
if rank == 0:
  
  model.plot_results() # plot results
  model.post_process()  # generate plots for the properties of MCMC sampling 

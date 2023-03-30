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
data = np.loadtxt("./sim_data.txt")
con = data[:126, :]
line= data[126:, :]

# make a data dict 
data_input = {"set1":[con, line]}
# if multiple datasets, e.g., 
#data_input = {"set1":[con1, line1], "set2":[con2, line2]}
# if a dataset has multiple lines, e.g.,
#data_input = {"set1":[con, line1, line2]}

#create a model
#there are two ways
#1) one way from the param file

#model = pymica.model(param_file="param/param_input")

#2) the ohter way is through the setup function

model = pymica.model()
# use Gaussians
model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[2, 2], max_num_saves=2000)

# or use tophats
#model.setup(data=data_input, type_tf='tophat', lag_limit=[0, 100], number_component=[2, 2], max_num_saves=2000)

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

#posterior run, only re-generate posterior samples, do not run MCMC
#model.post_run()

#do decomposition for the cases of multiple components 
#model.decompose()

# plot results
if rank == 0:
  
  model.plot_results() # plot results
  model.post_process()  # generate plots for the properties of MCMC sampling 

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
  timelag = model.get_posterior_sample_timelag(set=10, line=0) 
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
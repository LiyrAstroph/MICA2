# 
# an exemplary script to show the usage of pymica
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
  fname = "./sim_data.txt"
  data = np.loadtxt(fname)
  fp = open(fname)
  line = fp.readline()
  line = fp.readline()
  lsp = line[1:].split(":")
  
  nc = int(lsp[0])
  con = data[:nc, :]
  line= data[nc:, :]

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
# use gamma
model.setup(data=data_input, type_tf='exp', lag_limit=[0, 50], number_component=[2, 2], max_num_saves=2000)
# for gamma tf, lag_limit mean the time shift of the gamma function
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
#            flag_con_sys_err=False, flag_line_sys_err=False,
#            type_lag_prior=0, lag_prior=[[0, 50]],
#            num_particles=1, thread_steps_factor=1, 
#            new_level_interval_factor=1, save_interval_factor=1,
#            lam=10, beta=100, ptol=0.1, 
#            max_num_levels=0)

#run mica
model.run()

#resume from last run
# model.restart()

#posterior run, only re-generate posterior samples, do not run MCMC
#model.post_run()

#do decomposition for the cases of multiple components 
model.decompose()

# plot results
if rank == 0:
  
  model.plot_results(doshow=True, resp_input="resp_input.txt") # plot results, doshow controls whether showing the results on screen
  model.plot_decomp(doshow=True, resp_input="resp_input.txt")  # when model.decompose() is called
  model.post_process()  # generate plots for the properties of MCMC sampling 

  # get the full sample 
  # sample is a list, each element contains an array of posterior samples
  # sample[0] is for the case of number_component[0]
  # sample[1] is for the case of number_component[1] 
  # ...
  #sample = model.get_posterior_sample()

  # get the posterior sample of time lags of the "line" in the dataset "set"
  # timelag is a list, each element contains an array of posterior samples
  # timelag[0] is for the case of number_component[0]
  # timelag[1] is for the case of number_component[1]
  # ...
  #timelag = model.get_posterior_sample_timelag(set=0, line=0) 
  #plt.plot(timelag[0][:, 0])
  #plt.show()

  # get the posterior sample of widths of the "line" in the dataset "set"
  # width is a list, each element contains an array of posterior samples
  # width[0] is for the case of number_component[0]
  # width[1] is for the case of number_component[1]
  # ...
  #width = model.get_posterior_sample_width(set=0, line=0) 
  #plt.plot(width[0][:, 0])
  #plt.show()

# 
# an examplary script to show the useage of pymica
# 
# pmap module: photometric reverberation mapping
#
# the data is extracted from Fausnaugh et al. 2018, ApJ, 854, 10
#                   see ref. Ma, Q., Wu, X.-B. et al. 2023, arXiv:2303.10362
#
# run this script with "mpiexec -n np python example.py".
# replace "np" with the number of cores you want to use.
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
  con = np.loadtxt("mcg08_g.txt")
  line= np.loadtxt("mcg08_r.txt")
  
  # make a data dict 
  data_input = {"set1":[con, line]}
else:
  data_input = None 

data_input = comm.bcast(data_input, root=0)

#create a model
#there are two ways
#1) one way from the param file

#model = pymica.model(param_file="param/param_input")

#2) the ohter way is through the setup function
#   using pmap module.

model = pymica.pmap()
model.setup(data=data_input, type_tf='gaussian', max_num_saves=5000, lag_prior=[[-5, 5],[0, 50]], \
            ratio_prior=[[0.05, 0.5]], flag_con_sys_err=True, flag_line_sys_err=True, width_limit=[0.01, 100])

# if using top-hats, set type_tf='tophat'

#run mica
model.run()

#posterior run, only re-generate posterior samples, do not run MCMC
#model.post_run()

#do decomposition of the light curve into separate light curves of the two responses.
#model.decompose()

# plot results
if rank == 0:
  
  model.plot_results(doshow=True) # plot results,doshow controls whether showing the results on screen
  model.post_process()  # generate plots for the properties of MCMC sampling 

  model.plot_decomp()

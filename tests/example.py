from mpi4py import MPI
import numpy as np
import pymica
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

data = np.loadtxt("data/sim_data.txt")

con = data[:126, :]
line= data[126:, :]

data_input = {"set1":[con, line]}


#model = pymica.model(param_file="param/param")

model = pymica.model()
model.setup(data=data_input, type_tf='gaussian', lag_limit=[0, 100], number_component=[1, 1], max_num_saves=200)

model.postrun()

#model.decompose()

if rank == 0:
  model.plot_results()
  model.postprocess()

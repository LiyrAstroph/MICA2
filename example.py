from mpi4py import MPI
import numpy as np
import pymica
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

model = pymica.Basic()

model.run()
# Python wrapper of MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
# 

from __future__ import absolute_import

from libc.string cimport *  

from mpi4py import MPI
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

import numpy as np
cimport numpy as np
np.import_array()

from . import utility as ut

import os

cdef class model:
  """
  Model class
  """
  cdef PARSET parset
  cdef int rank, size

  def __cinit__(self, param_file=None):
    """
    initialise with a parameter file.
    """
    
    set_mpi()

    self.size  = MPI.COMM_WORLD.Get_size()
    self.rank  = MPI.COMM_WORLD.Get_rank()

    if param_file != None:
      if isinstance(param_file, str):
        strcpy(self.parset.param_file, param_file.encode("UTF-8"))
        set_param_file(self.parset.param_file)
        read_parset()
      else:
        raise ValueError("param_file should be a string!")
    
    if self.rank == 0:
      # check data folder
      if not os.path.exists("./data"):
        os.mkdir("./data/")
      # check param folder
      if not os.path.exists("./param"):
        os.mkdir("./param/")

    return
  
  def __cdealloc__(self):
    end_run()
    return
  

  def setup(self, data_file=None, data=None,
                  type_tf='gaussian', max_num_saves=2000, 
                  flag_uniform_var_params=False, flag_uniform_tranfuns=False,
                  flag_trend=0, flag_lag_posivity=False,
                  lag_limit=[0, 100], number_component=[1, 1],
                  flag_con_sys_err=False, flag_line_sys_err=False,
                  type_lag_prior=0):     
    """
    setup parameters
    """

    self._set_default_parset()

    # data file
    if data_file != None:
      if isinstance(data_file, str):
        strcpy(self.parset.data_file, data_file.encode("UTF-8"))
      else:
        raise ValueError("data_file is unrecognized!")
    elif data != None:
      self.create_data_file(data=data)
    else:
      raise ValueError("either data or data_file should not be None!")

    # transfer function type
    if type_tf == 'gaussian':
      self.parset.type_tf = 0
    elif type_tf == 'tophat':
      self.parset.type_tf = 1
    else:
      raise ValueError("type_tf is unrecognized!")
    
    # maximum number of saves
    self.parset.max_num_saves = max_num_saves

    # uniform variability parameters of mulitple datasets
    if flag_uniform_var_params == False:
      self.parset.flag_uniform_var_params = 0
    elif flag_uniform_var_params == True:
      self.parset.flag_uniform_var_params = 1
    else:
      raise ValueError("flag_uniform_var_params is unrecognized!")
    
    # uniform transfer function parameters of mulitple datasets
    if flag_uniform_tranfuns == False:
      self.parset.flag_uniform_tranfuns = 0
    elif flag_uniform_var_params == True:
      self.parset.flag_uniform_tranfuns = 1
    else:
      raise ValueError("flag_uniform_tranfuns is unrecognized!")
    
    # long-term trend
    self.parset.flag_trend = flag_trend 

    # lag posivitiy
    if flag_lag_posivity == False:
      self.parset.flag_lag_posivity = 0
    elif flag_lag_posivity == True:
      self.parset.flag_lag_posivity = 1
    else:
      raise ValueError("flag_lag_posivity is unrecognized!")

    # lag limit
    self.parset.lag_limit_low = lag_limit[0]
    self.parset.lag_limit_upper = lag_limit[1]

    # number of component
    self.parset.num_gaussian_low = number_component[0]
    self.parset.num_gaussian_upper = number_component[1]
    self.parset.num_gaussian_diff = self.parset.num_gaussian_upper-self.parset.num_gaussian_low + 1

    # continuum systematic error
    if flag_con_sys_err == False:
      self.parset.flag_con_sys_err = 0
    elif flag_con_sys_err == True:
      self.parset.flag_con_sys_err = 1
    else:
      raise ValueError("flag_lag_con_sys_err is unrecognized!")
    
    # line systematic error
    if flag_line_sys_err == False:
      self.parset.flag_line_sys_err = 0
    elif flag_line_sys_err == True:
      self.parset.flag_line_sys_err = 1
    else:
      raise ValueError("flag_lag_line_sys_err is unrecognized!")

    
    self.parset.type_lag_prior = type_lag_prior

    self.print_parset()

    # finally set parameters in C.
    set_parset(&self.parset)

    read_data()
    init()
    return
  
  def _set_default_parset(self):
    """
    set the default parameters
    """
    strcpy(self.parset.file_dir, './'.encode("UTF-8"))
    self.parset.num_gaussian_low = 1
    self.parset.num_gaussian_upper = 1
    self.parset.flag_uniform_tranfuns = 0
    self.parset.flag_uniform_var_params = 0
    self.parset.lag_limit_low = 0.0
    self.parset.lag_limit_upper = -1.0
    self.parset.type_lag_prior = 1
    self.parset.flag_trend = 0
    self.parset.type_tf = 0
    self.parset.flag_lag_posivity = 0
    self.parset.num_gaussian_low = 1
    self.parset.num_gaussian_upper = 1
    # cdnest options
    self.parset.num_particles = 2
    self.parset.max_num_saves = 2000
    self.parset.new_level_interval_factor = 2
    self.parset.save_interval_factor = self.parset.new_level_interval_factor
    self.parset.thread_steps_factor = 2
    self.parset.lam = 10.0
    self.parset.beta = 100.0
    self.parset.max_num_levels = 0
    self.parset.max_ptol = 0.1

    return
  
  def create_data_file(self, data=None):
    """
    create data file data_input.txt.
    """

    strcpy(self.parset.data_file, "data/data_input.txt".encode("UTF-8"))
    
    if self.rank != 0:
      return

    if data != None:
      if not isinstance(data, dict):
        raise ValueError("data is not a dict!")
      else:
        # check each datasets
        for key in data.keys():
          if len(data[key])<2:
            raise ValueError("dataset '%s' of data is not correct!"%key)

        fp = open(self.parset.data_file, "w")
        # write headers
        fp.write("# %d\n"%len(data))
        for key in data.keys():
          fp.write("# %d"%data[key][0].shape[0])
          for i in range(1, len(data[key])):
            fp.write(":%d"%data[key][i].shape[0])
          fp.write("\n")
        
        # write data
        for key in data.keys():
          for i in range(len(data[key])):
            np.savetxt(fp, data[key][i], fmt="%f")
            fp.write("\n")

        fp.close()

    return
  
  def print_parset(self):
    if self.rank == 0:
      fp = open("param/param_input", "w")
      fp.write("# lines starting with '#' are regarded as comments and are neglected\n")
      fp.write("# if want to turn on the line, remove the beginning '#'\n\n")
      fp.write("#==============================================================\n")

      fp.write("{:30}{}\n".format("FileDir", self.parset.file_dir.decode("UTF-8")))
      fp.write("{:30}{}\n".format("DataFile", self.parset.data_file.decode("UTF-8")))
      fp.write("{:30}{}\n".format("TypeTF", self.parset.type_tf))
      fp.write("{:30}{}\n".format("MaxNumberSaves", self.parset.max_num_saves))
      fp.write("{:30}{}\n".format("FlagUniformVarParams", self.parset.flag_uniform_var_params))
      fp.write("{:30}{}\n".format("FlagUniformTranFuns", self.parset.flag_uniform_tranfuns))
      fp.write("{:30}{}\n".format("FlagLongtermTrend", self.parset.flag_trend))
      fp.write("{:30}{}\n".format("LagLimitLow", self.parset.lag_limit_low))
      fp.write("{:30}{}\n".format("LagLimitUpp", self.parset.lag_limit_upper))
      fp.write("{:30}{}\n".format("FlagLagPositivity", self.parset.flag_lag_posivity))
      fp.write("{:30}{}\n".format("NumCompLow", self.parset.num_gaussian_low))
      fp.write("{:30}{}\n".format("NumCompUpp", self.parset.num_gaussian_upper))
      fp.write("{:30}{}\n".format("FlagConSysErr", self.parset.flag_con_sys_err))
      fp.write("{:30}{}\n".format("FlagLineSysErr", self.parset.flag_line_sys_err))
      fp.write("{:30}{}\n".format("TypeLagPrior", self.parset.type_lag_prior))
      fp.write("#==============================================================\n")
      fp.write("# options for cdnest sampling\n")
      fp.write("# use the default values or do not turn thme on IF NOT familiar with them\n\n")

      fp.write("{:30}{:5}{}\n".format("# PTol", self.parset.max_ptol,"# likelihood tolerance in loge"))
      fp.write("{:30}{:5}{}\n".format("# NumberParticles", self.parset.num_particles, "# number of particles"))
      fp.write("{:30}{:5}{}\n".format("# NewLevelIntervalFactor", self.parset.new_level_interval_factor, "# new level interval"))
      fp.write("{:30}{:5}{}\n".format("# SaveIntervalFactor",self.parset.save_interval_factor,"# particular saving interval"))
      fp.write("{:30}{:5}{}\n".format("# ThreadStepsFactor",self.parset.thread_steps_factor,"# thread steps before communications between cores"))
      fp.write("{:30}{:5}{}\n".format("# MaxNumberLevels",self.parset.max_num_levels,"# maximum number of levels; unlimited for 0"))
      fp.write("{:30}{:5}{}\n".format("# BacktrackingLength",self.parset.lam,"# backforward tracking length (lambda)"))
      fp.write("{:30}{:5}{}\n".format("# StrengthEqualPush",self.parset.beta, "# strength to force equal push (beta)"))
      fp.write("#==============================================================\n")
      fp.close()
    return

  def run(self):
    """
    run mica
    """
    mc_con()
    mc_line()
    
    return
  
  def postrun(self):
    """
    do posterior running
    """
    set_argv(1, 0)
    mc_line()
    return
  
  def plot_results(self):

    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      ut.plot_results(self.parset.file_dir.decode("UTF-8"), self.parset.data_file.decode("UTF-8"), i, \
                      self.parset.lag_limit_low, self.parset.lag_limit_upper, \
                      self.parset.flag_uniform_var_params, self.parset.flag_uniform_tranfuns, \
                      self.parset.flag_trend, self.parset.type_tf, None)
    return
  
  def postprocess(self, temperature=1.0):

    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      ut.postprocess(self.parset.file_dir.decode("UTF-8"), i, temperature)

    return
  
  def decompose(self):
    """
    do decomposition
    """
    set_argv(1, 1)
    mc_line()
    return
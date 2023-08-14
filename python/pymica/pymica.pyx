# Python wrapper of MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
# 

from __future__ import absolute_import

from libc.string cimport *  
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc

from mpi4py import MPI
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

import numpy as np
cimport numpy as np
np.import_array()

from . import utility as ut

import os

__all__ = ["gmodel", "pmap", "dmap"]

cdef class basis:

  cdef PARSET parset
  cdef int rank, size
  cdef int nset, num_param_var 
  cdef list nlset

  def __cinit__(self):
    
    set_mpi()  # setup MPI enviroment in C
    self._set_default_parset()

    self.size  = MPI.COMM_WORLD.Get_size()
    self.rank  = MPI.COMM_WORLD.Get_rank()

    self.nset = 0
    self.num_param_var = 0

    if self.rank == 0:
      # check data folder
      if not os.path.exists("./data"):
        os.mkdir("./data/")
      # check param folder
      if not os.path.exists("./param"):
        os.mkdir("./param/")

    return

  def __cdealloc__(self):
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
    self.parset.width_limit_low_isset = 0
    self.parset.width_limit_upper_isset = 0
    self.parset.type_lag_prior = 1
    self.parset.flag_trend = 0
    self.parset.type_tf = 0
    self.parset.flag_lag_posivity = 0
    self.parset.num_gaussian_low = 1
    self.parset.num_gaussian_upper = 1
    # cdnest options
    self.parset.num_particles = 1
    self.parset.max_num_saves = 2000
    self.parset.new_level_interval_factor = 1
    self.parset.save_interval_factor = self.parset.new_level_interval_factor
    self.parset.thread_steps_factor = 1
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
      fp.write("{:30}{}\n".format("TypeModel", self.parset.model))
      fp.write("{:30}{}\n".format("TypeTF", self.parset.type_tf))
      fp.write("{:30}{}\n".format("MaxNumberSaves", self.parset.max_num_saves))
      fp.write("{:30}{}\n".format("FlagUniformVarParams", self.parset.flag_uniform_var_params))
      fp.write("{:30}{}\n".format("FlagUniformTranFuns", self.parset.flag_uniform_tranfuns))
      fp.write("{:30}{}\n".format("FlagLongtermTrend", self.parset.flag_trend))
      if self.parset.model == 0 or self.parset.model == 2:
        fp.write("{:30}{}\n".format("LagLimitLow", self.parset.lag_limit_low))
        fp.write("{:30}{}\n".format("LagLimitUpp", self.parset.lag_limit_upper))
      if self.parset.width_limit_low_isset == 1:
        fp.write("{:30}{}\n".format("WidthLimitLow", self.parset.width_limit_low))
      if self.parset.width_limit_upper_isset == 1:
        fp.write("{:30}{}\n".format("WidthLimitUpp", self.parset.width_limit_upper))  
      fp.write("{:30}{}\n".format("FlagLagPositivity", self.parset.flag_lag_posivity))
      fp.write("{:30}{}\n".format("NumCompLow", self.parset.num_gaussian_low))
      fp.write("{:30}{}\n".format("NumCompUpp", self.parset.num_gaussian_upper))
      fp.write("{:30}{}\n".format("FlagConSysErr", self.parset.flag_con_sys_err))
      fp.write("{:30}{}\n".format("FlagLineSysErr", self.parset.flag_line_sys_err))
      fp.write("{:30}{}\n".format("TypeLagPrior", self.parset.type_lag_prior))
      if self.parset.type_lag_prior == 4:
        fp.write("{:30}{}\n".format("StrLagPrior", self.parset.str_lag_prior.decode("UTF-8")))
      if strlen(self.parset.str_ratio_prior)>0:
        fp.write("{:30}{}\n".format("StrRatioPrior", self.parset.str_ratio_prior.decode("UTF-8")))
      fp.write("#==============================================================\n")
      fp.write("# options for cdnest sampling\n")
      fp.write("# use the default values or do not turn them on IF NOT familiar with them\n\n")

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
  
  def plot_results(self):

    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      ut.plot_results(self.parset.file_dir.decode("UTF-8"), self.parset.data_file.decode("UTF-8"), i, \
                      self.parset.lag_limit_low, self.parset.lag_limit_upper, \
                      self.parset.flag_uniform_var_params, self.parset.flag_uniform_tranfuns, \
                      self.parset.flag_trend, self.parset.type_tf, self.parset.model, None)
    return
  
  def plot_decomp(self):
    """
    plot line decomposition
    """
    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      ut.plot_line_decomp(self.parset.file_dir.decode("UTF-8"), self.parset.data_file.decode("UTF-8"), i, \
                          self.parset.type_tf, self.parset.model)
    
    return
  
  def post_process(self, temperature=1.0):

    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      ut.postprocess(self.parset.file_dir.decode("UTF-8"), i, temperature)

    return
  
  def get_posterior_sample_width(self, set=0, line=0):
    """
    get the posterior sample of widths of the line in set  
    """  

    if set >= self.nset:
      raise ValueError("there are only %d datasets!\n set exceeds this limit!"%self.nset)
    elif set < 0:
      raise ValueError("set should be non-negative!")
    
    if line >= self.nlset[set]:
      raise ValueError("there are only %d lines!\n line exceeds this limit!"%self.nset)
    elif line < 0:
      raise ValueError("line should be non-negative!")

    width = []
    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      sample = np.loadtxt(self.parset.file_dir.decode("UTF-8")+"/data/posterior_sample1d.txt_%d"%i)
      if self.parset.flag_uniform_tranfuns == 0:
        idx_line = self.num_param_var
        for j in range(0, set-1):
          idx_line += (1+(i*3))*self.nlset[j]

        idx_line += (1+(i*3))*line
        width.append(sample[:, idx_line+3:idx_line+3+i*3:3])
      else:
        idx_line = self.num_param_var
        idx_line += (1+(i*3))*line
        width.append(sample[:, idx_line+3:idx_line+3+i*3:3])

    return width
  
  def get_posterior_sample(self):
    """
    get posterior sample
    """
    sample=[]
    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      sample.append(np.loadtxt(self.parset.file_dir.decode("UTF-8")+"/data/posterior_sample1d.txt_%d"%i))
    
    return sample 
  
  def get_posterior_sample_timelag(self, set=0, line=0):
    """
    get the posterior sample of time lags of the line in set  
    """    
    if set >= self.nset:
      raise ValueError("there are only %d datasets!\n set exceeds this limit!"%self.nset)
    elif set < 0:
      raise ValueError("set should be non-negative!")
    
    if line >= self.nlset[set]:
      raise ValueError("there are only %d lines!\n line exceeds this limit!"%self.nset)
    elif line < 0:
      raise ValueError("line should be non-negative!")

    timelag = []
    for i in range(self.parset.num_gaussian_low, self.parset.num_gaussian_upper+1, 1):
      sample = np.loadtxt(self.parset.file_dir.decode("UTF-8")+"/data/posterior_sample1d.txt_%d"%i)
      if self.parset.flag_uniform_tranfuns == 0:
        idx_line = self.num_param_var
        for j in range(0, set-1):
          idx_line += (1+(i*3))*self.nlset[j]

        idx_line += (1+(i*3))*line
        timelag.append(sample[:, idx_line+2:idx_line+2+i*3:3])
      else:
        idx_line = self.num_param_var
        idx_line += (1+(i*3))*line
        timelag.append(sample[:, idx_line+2:idx_line+2+i*3:3])

    return timelag
  
  def _get_data_dimension(self):
    """
    get data dimensions
    """
    if self.rank == 0:
      fp = open(self.parset.data_file.decode("UTF-8"), "r")
      line = fp.readline()
      self.nset = int(line[1:])
      self.nlset = []
      for i in range(self.nset):
        line = fp.readline()
        ls = line[1:].split(":")
        self.nlset.append(len(ls)-1)
      
      fp.close()
    
    self.nset = MPI.COMM_WORLD.bcast(self.nset, root=0)
    self.nlset = MPI.COMM_WORLD.bcast(self.nlset, root=0)
    
    return
  
  def setup(self, data_file=None, data=None,
                  type_tf='gaussian', max_num_saves=2000,
                  flag_trend=0, flag_lag_posivity=False,
                  flag_con_sys_err=False, flag_line_sys_err=False,
                  # follows cdnest parameters
                  num_particles=1, thread_steps_factor=1, 
                  new_level_interval_factor=1, save_interval_factor=1,
                  lam=10, beta=100, ptol=0.1, 
                  max_num_levels=0):
    """
    setup parameters
    """
    # data file
    if data_file != None:
      if isinstance(data_file, str):
        strcpy(self.parset.data_file, data_file.encode("UTF-8"))
        self._get_data_dimension()
      else:
        raise ValueError("data_file is unrecognized!")
    elif data != None:
      self.create_data_file(data=data)
      # set data dimensions
      self.nset = len(data)
      self.nlset = []
      for key in data.keys():
        self.nlset.append(len(data[key])-1)

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

    # long-term trend
    self.parset.flag_trend = flag_trend 

    # lag posivitiy
    if flag_lag_posivity == False:
      self.parset.flag_lag_posivity = 0
    elif flag_lag_posivity == True:
      self.parset.flag_lag_posivity = 1
    else:
      raise ValueError("flag_lag_posivity is unrecognized!")
    
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
    
    self.parset.num_particles = num_particles
    self.parset.thread_steps_factor = thread_steps_factor
    self.parset.new_level_interval_factor = new_level_interval_factor
    self.parset.save_interval_factor = save_interval_factor
    self.parset.max_num_levels = max_num_levels
    self.parset.lam = lam
    self.parset.beta = beta
    self.parset.max_ptol = ptol

    return


#==========================================================================
# general model 
cdef class gmodel(basis):
  """
  Model class
  """
  def __cinit__(self, param_file=None):
    """
    initialise with a parameter file.
    """
    if param_file != None:
      if isinstance(param_file, str):
        strcpy(self.parset.param_file, param_file.encode("UTF-8"))
        set_param_file(self.parset.param_file)
        read_parset()
        get_parset(&self.parset)  # get parset from C

        self._get_data_dimension()
        # uniform variability parameters of mulitple datasets
        if self.parset.flag_uniform_var_params == 0:
          self.num_param_var = 3*self.nset
        elif self.parset.flag_uniform_var_params == 1:
          self.num_param_var = 3

      else:
        raise ValueError("param_file should be a string!")

    return
  
  def __cdealloc__(self):
    return
  

  def setup(self, data_file=None, data=None,
                  type_tf='gaussian', max_num_saves=2000, 
                  flag_uniform_var_params=False, flag_uniform_tranfuns=False,
                  flag_trend=0, flag_lag_posivity=False,
                  lag_limit=[0, 100], number_component=[1, 1],
                  width_limit=None,
                  flag_con_sys_err=False, flag_line_sys_err=False,
                  type_lag_prior=0, lag_prior=None,
                  # follows cdnest parameters
                  num_particles=1, thread_steps_factor=1, 
                  new_level_interval_factor=1, save_interval_factor=1,
                  lam=10, beta=100, ptol=0.1, 
                  max_num_levels=0):     
    """
    setup parameters
    """
    
    basis.setup(self, data_file, data, type_tf, max_num_saves, flag_trend, flag_lag_posivity, \
                      flag_con_sys_err, flag_line_sys_err,
                      # follows cdnest parameters
                      num_particles = num_particles,
                      thread_steps_factor = thread_steps_factor, 
                      new_level_interval_factor = new_level_interval_factor,
                      save_interval_factor = save_interval_factor,
                      lam = lam,
                      beta = beta, 
                      ptol = ptol, 
                      max_num_levels = max_num_levels)

    # uniform variability parameters of mulitple datasets
    if flag_uniform_var_params == False:
      self.parset.flag_uniform_var_params = 0
      self.num_param_var = 3*self.nset
    elif flag_uniform_var_params == True:
      self.parset.flag_uniform_var_params = 1
      self.num_param_var = 3
    else:
      raise ValueError("flag_uniform_var_params is unrecognized!")
    
    # uniform transfer function parameters of mulitple datasets
    if flag_uniform_tranfuns == False:
      self.parset.flag_uniform_tranfuns = 0
    elif flag_uniform_tranfuns == True:
      self.parset.flag_uniform_tranfuns = 1
    else:
      raise ValueError("flag_uniform_tranfuns is unrecognized!")

    # lag limit
    self.parset.lag_limit_low = lag_limit[0]
    self.parset.lag_limit_upper = lag_limit[1]

    # width limit
    if width_limit != None:
      self.parset.width_limit_low_isset = 1
      self.parset.width_limit_low = width_limit[0]
      self.parset.width_limit_upper_isset = 1
      self.parset.width_limit_upper = width_limit[1]

    # number of component
    if isinstance(number_component, int):
      self.parset.num_gaussian_low = number_component
      self.parset.num_gaussian_upper = number_component
    elif isinstance(number_component, list):
      self.parset.num_gaussian_low = number_component[0]
      self.parset.num_gaussian_upper = number_component[1]
    self.parset.num_gaussian_diff = self.parset.num_gaussian_upper-self.parset.num_gaussian_low + 1
    
    self.parset.type_lag_prior = type_lag_prior

    # if lag_prior is input
    if lag_prior != None:
      if self.parset.num_gaussian_upper > 1:
        self.parset.type_lag_prior=4
        # write string of lag prior
        sstr = "["
        sstr += "%f:%f"%(lag_prior[0][0], lag_prior[0][1])
        for i in range(1, self.parset.num_gaussian_upper):
          sstr += ":%f:%f"%(lag_prior[i][0], lag_prior[i][1])
        
        sstr += "]"
        strcpy(self.parset.str_lag_prior, sstr.encode("UTF-8"))
      else:
        self.parset.type_lag_prior=0

    self.print_parset()

    # finally set parameters in C.
    set_parset(&self.parset)
    return
  
  def run(self):
    """
    run mica
    """
    read_data()
    init()
    mc_con()
    mc_line()
    end_run()
    return
  
  def post_run(self):
    """
    do posterior running
    """
    set_argv(1, 0, 0) # postprocess, decompose, restart
    read_data()
    init()
    mc_line()
    end_run()
    return
  
  def decompose(self):
    """
    do decomposition
    """
    set_argv(1, 1, 0) # postprocess, decompose, restart
    read_data()
    init()
    mc_line()
    end_run()
    return
  
  def restart(self): 
    """
    resume from a last run
    """
    set_argv(0, 0, 1)
    read_data()
    init()
    mc_line()
    end_run()
    return


#==========================================================================
# photometric reverberation mapping
#
cdef class pmap(basis):

  def __cinit__(self):
    return 
  
  def __cdealloc__(self):
    return
  
  def setup(self, data_file=None, data=None,
                  type_tf='gaussian', max_num_saves=2000,
                  flag_trend=0, flag_lag_posivity=False,
                  flag_con_sys_err=False, flag_line_sys_err=False,
                  lag_prior=None, ratio_prior=None,
                  width_limit=None,
                  # follows cdnest parameters
                  num_particles=1, thread_steps_factor=1, 
                  new_level_interval_factor=1, save_interval_factor=1,
                  lam=10, beta=100, ptol=0.1, 
                  max_num_levels=0):

    basis.setup(self, data_file, data, type_tf, max_num_saves, flag_trend, flag_lag_posivity, \
                      flag_con_sys_err, flag_line_sys_err,
                      # follows cdnest parameters
                      num_particles = num_particles,
                      thread_steps_factor = thread_steps_factor, 
                      new_level_interval_factor = new_level_interval_factor,
                      save_interval_factor = save_interval_factor,
                      lam = lam,
                      beta = beta, 
                      ptol = ptol, 
                      max_num_levels = max_num_levels)
    
    self.parset.num_gaussian_low = 2
    self.parset.num_gaussian_upper = 2 
    self.parset.num_gaussian_diff = self.parset.num_gaussian_upper - self.parset.num_gaussian_low + 1
    
    self.parset.model=1
    self.parset.type_lag_prior=4
    if lag_prior == None:
      raise ValueError("Please input lag_prior.\ne.g., lag_prior=[[-5, 5], [5, 500]]\n")
    else:
      # write string of lag prior
      sstr = "["
      sstr += "%f:%f"%(lag_prior[0][0], lag_prior[0][1])
      for i in range(1, self.parset.num_gaussian_upper):
        sstr += ":%f:%f"%(lag_prior[i][0], lag_prior[i][1])
      
      sstr += "]"
      strcpy(self.parset.str_lag_prior, sstr.encode("UTF-8"))
   
    # if lag_prior is input
    if ratio_prior != None:
      
      if not isinstance(ratio_prior[0], list):
        ratio_prior = [ratio_prior]

      # write string of lag prior
      sstr = "["
      sstr += "%f:%f"%(ratio_prior[0][0], ratio_prior[0][1])
      for i in range(1, len(ratio_prior)):
        sstr += ":%f:%f"%(ratio_prior[i][0], ratio_prior[i][1])
      
      sstr += "]"
      strcpy(self.parset.str_ratio_prior, sstr.encode("UTF-8"))

    # width limit
    if width_limit != None:
      self.parset.width_limit_low_isset = 1
      self.parset.width_limit_low = width_limit[0]
      self.parset.width_limit_upper_isset = 1
      self.parset.width_limit_upper = width_limit[1]

    self.print_parset()

    # finally set parameters in C.
    set_parset(&self.parset)
    
    return
  
  def run(self):
    """
    run mica
    """
    read_data()
    init()
    mc_con()
    mc_pmap()
    end_run()
    return
  
  def post_run(self):
    """
    do posterior running
    """
    set_argv(1, 0, 0) # postprocess, decompose, restart
    read_data()
    init()
    mc_pmap()
    end_run()
    return
  
  def decompose(self):
    """
    do decomposition
    """
    set_argv(1, 1, 0)  # postprocess, decompose, restart
    read_data()
    init()
    mc_pmap()
    end_run()
    return
  
  def restart(self): 
    """
    resume from a last run
    """
    set_argv(0, 0, 1)
    read_data()
    init()
    mc_pmap()
    end_run()
    return


#==========================================================================
# disk reverberation mapping 
# the continuum data is empty
cdef class dmap(basis):
  """
  Model class
  """
  def __cinit__(self, param_file=None):
    """
    initialise with a parameter file.
    """
    if param_file != None:
      if isinstance(param_file, str):
        strcpy(self.parset.param_file, param_file.encode("UTF-8"))
        set_param_file(self.parset.param_file)
        read_parset()
        get_parset(&self.parset)  # get parset from C

        self._get_data_dimension()
        # uniform variability parameters of mulitple datasets
        if self.parset.flag_uniform_var_params == 0:
          self.num_param_var = 3*self.nset
        elif self.parset.flag_uniform_var_params == 1:
          self.num_param_var = 3

      else:
        raise ValueError("param_file should be a string!")

    return
  
  def __cdealloc__(self):
    return
  

  def setup(self, data_file=None, data=None,
                  type_tf='gaussian', max_num_saves=2000, 
                  flag_uniform_var_params=False, flag_uniform_tranfuns=False,
                  flag_trend=0, flag_lag_posivity=False,
                  lag_limit=[0, 100], number_component=[1, 1],
                  width_limit=None,
                  flag_con_sys_err=False, flag_line_sys_err=False,
                  type_lag_prior=0, lag_prior=None,
                  # follows cdnest parameters
                  num_particles=1, thread_steps_factor=1, 
                  new_level_interval_factor=1, save_interval_factor=1,
                  lam=10, beta=100, ptol=0.1, 
                  max_num_levels=0):     
    """
    setup parameters
    """
    
    basis.setup(self, data_file, data, type_tf, max_num_saves, flag_trend, flag_lag_posivity, \
                      flag_con_sys_err, flag_line_sys_err,
                      # follows cdnest parameters
                      num_particles = num_particles,
                      thread_steps_factor = thread_steps_factor, 
                      new_level_interval_factor = new_level_interval_factor,
                      save_interval_factor = save_interval_factor,
                      lam = lam,
                      beta = beta, 
                      ptol = ptol, 
                      max_num_levels = max_num_levels)
   
    self.parset.model = 2

    # uniform variability parameters of mulitple datasets
    if flag_uniform_var_params == False:
      self.parset.flag_uniform_var_params = 0
      self.num_param_var = 3*self.nset
    elif flag_uniform_var_params == True:
      self.parset.flag_uniform_var_params = 1
      self.num_param_var = 3
    else:
      raise ValueError("flag_uniform_var_params is unrecognized!")
    
    # uniform transfer function parameters of mulitple datasets
    if flag_uniform_tranfuns == False:
      self.parset.flag_uniform_tranfuns = 0
    elif flag_uniform_tranfuns == True:
      self.parset.flag_uniform_tranfuns = 1
    else:
      raise ValueError("flag_uniform_tranfuns is unrecognized!")

    # lag limit
    self.parset.lag_limit_low = lag_limit[0]
    self.parset.lag_limit_upper = lag_limit[1]

    # width limit
    if width_limit != None:
      self.parset.width_limit_low_isset = 1
      self.parset.width_limit_low = width_limit[0]
      self.parset.width_limit_upper_isset = 1
      self.parset.width_limit_upper = width_limit[1]

    # number of component
    if isinstance(number_component, int):
      self.parset.num_gaussian_low = number_component
      self.parset.num_gaussian_upper = number_component
    elif isinstance(number_component, list):
      self.parset.num_gaussian_low = number_component[0]
      self.parset.num_gaussian_upper = number_component[1]
    self.parset.num_gaussian_diff = self.parset.num_gaussian_upper-self.parset.num_gaussian_low + 1
    
    self.parset.type_lag_prior = type_lag_prior

    # if lag_prior is input
    if lag_prior != None:
      if self.parset.num_gaussian_upper > 1:
        self.parset.type_lag_prior=4
        # write string of lag prior
        sstr = "["
        sstr += "%f:%f"%(lag_prior[0][0], lag_prior[0][1])
        for i in range(1, self.parset.num_gaussian_upper):
          sstr += ":%f:%f"%(lag_prior[i][0], lag_prior[i][1])
        
        sstr += "]"
        strcpy(self.parset.str_lag_prior, sstr.encode("UTF-8"))
      else:
        self.parset.type_lag_prior=0

    self.print_parset()

    # finally set parameters in C.
    set_parset(&self.parset)
    return
  
  def run(self):
    """
    run mica
    """
    read_data()
    init()
    mc_dmap()
    end_run()
    return
  
  def post_run(self):
    """
    do posterior running
    """
    set_argv(1, 0, 0) # postprocess, decompose, restart
    read_data()
    init()
    mc_dmap()
    end_run()
    return
  
  def decompose(self):
    """
    do decomposition
    """
    set_argv(1, 1, 0) # postprocess, decompose, restart
    read_data()
    init()
    mc_dmap()
    end_run()
    return
  
  def restart(self): 
    """
    resume from a last run
    """
    set_argv(0, 0, 1)
    read_data()
    init()
    mc_dmap()
    end_run()
    return

#!/usr/bin/python
#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=True, profile=False

#
# Python wrapper of MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
# 

from libc.stdio cimport * 

cdef extern from "allvars.h":
  ctypedef struct PARSET:
    char param_file[256]
    char file_dir[256]
    char data_file[256]
    char str_lag_prior[256]
    char str_ratio_prior[256]
  
    int flag_uniform_var_params
    int flag_uniform_tranfuns
    int flag_trend
    int flag_con_sys_err, flag_line_sys_err
    int type_lag_prior
    int type_tf
    int model
    int flag_lag_posivity
  
    int num_gaussian_low, num_gaussian_upper, num_gaussian_diff
  
    double lag_limit_low, lag_limit_upper
    int width_limit_low_isset, width_limit_upper_isset
    double width_limit_low, width_limit_upper
    
    # otpions for cdnest 
    int max_num_saves
    int num_particles
    int max_num_levels
    double lam, beta, max_ptol
    double thread_steps_factor, new_level_interval_factor, save_interval_factor
  
cdef extern from "proto.h":
  void set_mpi()
  void set_param_file(char *fname)
  void set_parset(PARSET *ps)
  void get_parset(PARSET *ps)
  void read_parset()
  void read_data()
  void init()
  void mc_con()
  void mc_line()
  void mc_pmap()
  void mc_vmap()
  void output_decompose_line()
  void begin_run()
  void end_run()
  void set_argv(int flag_pp, int flag_de, int flag_re)

/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef MICA_ALLVARS_H
#define MICA_ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define PI            (M_PI)
#define MICA_MAX_STR_LENGTH  (256)

/* variables for MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];


/* for dnest */
extern int which_level_update;
extern double *limits;  // external from dnest
extern char dnest_options_file[MICA_MAX_STR_LENGTH], 
            dnest_options_file_con[MICA_MAX_STR_LENGTH],
            dnest_options_file_line[MICA_MAX_STR_LENGTH];
extern char postfix[MICA_MAX_STR_LENGTH];

/* con */
extern int num_params, num_params_var, num_params_line;
extern double **par_range_model;
extern int *par_fix;
extern double *par_fix_val;
extern double var_range_model[7][2], line_range_model[4][2];
extern double *PNmat, *PSmat, *PCmat, *IPCmat, *USmat, *USmatT, *ASmat, *Tmat1, *Tmat2;
extern double *workspace;
extern int ncon_max, nall_max, nlset_max;
extern int nq; /* order of long-term trend */
extern double *var_param, *var_param_std;
extern double tspan_max, tcadence_con_min, tcadence_line_min, tcadence_min;


/* line */
extern int **idx_line_pm, *idx_con_pm;
extern int num_gaussian;
extern int nline_max;
extern int type_lag_prior_pr;
extern double *Smat_lc, *Smat_ll;
extern double *logz_arr; /* evidence */

// error exit
extern char str_error_exit[200];
// mathematic functions
extern int *workspace_ipiv;
extern double *workspace_inv;


extern int flag_postprc, flag_end, flag_restart, flag_decomp;

typedef struct 
{
  char param_file[MICA_MAX_STR_LENGTH];
  char file_dir[MICA_MAX_STR_LENGTH];
  char data_file[MICA_MAX_STR_LENGTH];

  int flag_uniform_var_params;
  int flag_uniform_tranfuns;
  int flag_trend;
  int flag_con_sys_err, flag_line_sys_err;
  int type_lag_prior;
  int type_tf;
  int flag_lag_posivity;

  int num_gaussian_low, num_gaussian_upper, num_gaussian_diff;

  double lag_limit_low, lag_limit_upper;
  
  /* otpions for cdnest */
  int max_num_saves;
  int num_particles;
  int max_num_levels;
  double lam, beta, max_ptol;
  double thread_steps_factor, new_level_interval_factor, save_interval_factor;
}PARSET;
extern PARSET parset;

typedef struct
{
  int id;
  void *addr;
  char tag[50];
  int isset;
}PARDICT;
extern PARDICT *pardict;
extern int num_pardict;

typedef struct
{
  int n;
  double scale, error_mean;
  double *t, *f, *fe;
}LC;

typedef struct /* each dataset contains one continuum and several line */
{
  int nlset;
  LC con;
  LC *line;
}DATASET;

extern DATASET *dataset; 
extern int nset;

extern LC *alldata;

extern int nscale, nrec_max; /* determine number of points of reconstruction */

// whether print warning
extern int flag_warning;
#endif
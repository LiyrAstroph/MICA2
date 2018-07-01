/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define PI            (M_PI)
#define MICA_MAX_STR_LENGTH  (100)

/* variables for MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];


/* for dnest */
extern int which_particle_update;
extern int which_level_update;
extern unsigned long long int which_mcmc_steps;//mcmc steps 
extern int *perturb_accept;
extern double *limits;  // external from dnest
extern char dnest_options_file[MICA_MAX_STR_LENGTH];

/* con */
extern int num_params, num_params_var, num_params_line;
extern double **par_range_model;
extern int *par_fix;
extern double *par_fix_val;
extern double var_range_model[7][2], line_range_model[4][2];
extern double *PNmat, *PSmat, *PCmat, *IPCmat, *USmat, *ASmat;
extern double *workspace;
extern int ncon_max, nall_max, nlset_max;
extern int nq; /* order of long-term trend */
extern double *var_param, *var_param_std;


/* line */
extern int **idx_line_pm;

// error exit
extern char str_error_exit[200];
// mathematic functions
extern int *workspace_ipiv;

typedef struct 
{
  char param_file[MICA_MAX_STR_LENGTH];
  char file_dir[MICA_MAX_STR_LENGTH];
  char data_file[MICA_MAX_STR_LENGTH];

  int num_particles;

  int flag_uniform_var_params;
  int flag_uniform_tranfuns;

}PARSET;
extern PARSET parset;

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

#endif
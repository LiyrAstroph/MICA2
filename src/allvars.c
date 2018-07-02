/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"

/* MPICH */
int thistask, totaltask, namelen;
int roottask = 0;
char proc_name[MPI_MAX_PROCESSOR_NAME];

/* dnest */
int which_particle_update;
int which_level_update;
unsigned long long int which_mcmc_steps;/* mcmc steps */ 
int *perturb_accept;
char dnest_options_file[MICA_MAX_STR_LENGTH];

PARSET parset;

DATASET *dataset;
int nset;

LC *alldata;

/* con */
int num_params, num_params_var, num_params_line;
double **par_range_model;
int *par_fix;
double *par_fix_val;
double var_range_model[7][2], line_range_model[4][2];;
double *PNmat, *PSmat, *PCmat, *IPCmat, *USmat, *ASmat, *Tmat1, *Tmat2;
double *workspace;
int ncon_max, nall_max, nlset_max;
int nq; /* order of long-term trend */
double *var_param, *var_param_std;

/* line */
int **idx_line_pm, *idx_con_pm; 

// mathematic functions
int *workspace_ipiv;

// error exit
char str_error_exit[200];

/* GSL */
gsl_integration_fixed_workspace * gsl_w;
gsl_integration_fixed_type * gsl_T;
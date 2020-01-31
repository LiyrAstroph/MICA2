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

#include "allvars.h"

/* MPICH */
int thistask, totaltask, namelen;
int roottask = 0;
char proc_name[MPI_MAX_PROCESSOR_NAME];

/* dnest */
int which_level_update;
char dnest_options_file[MICA_MAX_STR_LENGTH];
char postfix[MICA_MAX_STR_LENGTH];

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
double *PNmat, *PSmat, *PCmat, *IPCmat, *USmat, *USmatT, *ASmat, *Tmat1, *Tmat2;
double *workspace;
int ncon_max, nall_max, nlset_max;
int nq; /* order of long-term trend */
double *var_param, *var_param_std;
double tspan_max, tcadence_con_min, tcadence_line_min, tcadence_min;

/* line */
int **idx_line_pm, *idx_con_pm; 
int num_gaussian;
int nline_max;
int type_lag_prior_pr;
double *Smat_lc, *Smat_ll;
double *logz_arr; /* evidence */

// mathematic functions
int *workspace_ipiv;
double *workspace_inv;

// error exit
char str_error_exit[200];

int flag_postprc;

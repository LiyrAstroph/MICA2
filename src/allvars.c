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
char dnest_options_file[MICA_MAX_STR_LENGTH], 
     dnest_options_file_con[MICA_MAX_STR_LENGTH],
     dnest_options_file_line[MICA_MAX_STR_LENGTH];
char postfix[MICA_MAX_STR_LENGTH];

PARSET parset;

DATASET *dataset;
int nset;

PARDICT *pardict;
int num_pardict;

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
double *lag_prior, *width_prior;

/* pmap */
void *model_trans;
double *resp_ratio_prior;

/* for seasonal gap */
double *gap_width, *gap_center;

/* for max and min flux of each dataset */
double **flux_minmax;

// mathematic functions
int *workspace_ipiv;
double *workspace_inv;

// error exit
char str_error_exit[200];

int flag_postprc, flag_end, flag_restart, flag_decomp, flag_examine, 
    flag_para_name, flag_load_prior;
char prior_file[MICA_MAX_STR_LENGTH];

int flag_warning = 0;

int nscale, nrec_max;

double width_factor, lag_factor;
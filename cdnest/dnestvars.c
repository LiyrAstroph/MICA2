/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "dnestvars.h"
 
/* output files */
FILE *fsample, *fsample_info;

/* random number generator */
const gsl_rng_type * dnest_gsl_T;
gsl_rng * dnest_gsl_r;

DNestOptions options;
char options_file[STR_MAX_LENGTH];

// sampler
bool save_to_disk;
double compression;
unsigned int regularisation;

void *particles;
int dnest_size_of_modeltype;
int particle_offset_size, particle_offset_double;

LikelihoodType *log_likelihoods;
unsigned int *level_assignments;

// number account of unaccepted times
unsigned int *account_unaccepts;

int size_levels, size_levels_combine;  
Level *levels;
Level *copies_of_levels, *levels_combine;
LikelihoodType *all_above;
unsigned int count_saves, num_saves, num_saves_restart;
int dnest_which_particle_update; // which particle to be updated
int dnest_which_level_update;    // which level to be updated;
unsigned long long int count_mcmc_steps;
LikelihoodType *above;
unsigned int size_above, size_all_above;

double post_logz;
int dnest_num_params;
char dnest_sample_postfix[STR_MAX_LENGTH], dnest_sample_tag[STR_MAX_LENGTH], dnest_sample_dir[STR_MAX_LENGTH];
double *dnest_param_range, *dnest_prior_info;
int *dnest_prior_type;
void *dnest_args;

double *limits, *copies_of_limits;

int dnest_thistask, dnest_totaltask;
int *dnest_perturb_accept;
int dnest_root;

int dnest_flag_restart=0, dnest_flag_postprc=0, dnest_flag_sample_info=0, dnest_flag_limits=0;
double dnest_post_temp=1.0;
char file_restart[STR_MAX_LENGTH], file_save_restart[STR_MAX_LENGTH];

//***********************************************
/*                  functions                  */
double mod(double y, double x);
void wrap_limit(double *x, double min, double max);
int mod_int(int y, int x);
int dnest_cmp(const void *pa, const void *pb);

void options_load(char *optfile, DNestOptions *opts);
void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, 
           double *param_range, int *prior_type, double *prior_info, 
           char *sample_dir, char *optfile, DNestOptions *opts, void *args);
void finalise();

void dnest_run();
void dnest_mcmc_run();
void update_particle(unsigned int which);
void update_level_assignment(unsigned int which);
double log_push(unsigned int which_level);
bool enough_levels(Level *l, int size_l);
void do_bookkeeping();
void save_levels();
void save_particle();
void save_limits();
void kill_lagging_particles();
void renormalise_visits();
void recalculate_log_X();
void dnest_postprocess(double temperature, char *optfile, DNestOptions *opts);
void postprocess(double temperature);
void initialize_output_file();
void close_output_file();
void dnest_save_restart();
void dnest_restart();
void dnest_restart_action(int iflag);
void dnest_accept_action();
void dnest_kill_action(int i, int i_copy);
void dnest_from_prior(void *model);
double dnest_perturb(void *model);
void dnest_print_particle(FILE *fp, const void *model);
void dnest_read_particle(FILE *fp, void *model);
/*=====================================================*/
// users responsible for following functions
void (*print_particle)(FILE *fp, const void *model);
void (*read_particle)(FILE *fp, void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
void (*restart_action)(int iflag);
void (*accept_action)();
void (*kill_action)(int i, int i_copy);
/*=====================================================*/
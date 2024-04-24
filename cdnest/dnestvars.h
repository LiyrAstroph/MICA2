/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */
#ifndef _DNESTVARS_H
#define _DNESTVARS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "dnest.h"

enum PRIOR_TYPE {UNIFORM=0, GAUSSIAN=1, LOG=2};

/* output files */
extern FILE *fsample, *fsample_info;

/* random number generator */
extern const gsl_rng_type * dnest_gsl_T;
extern gsl_rng * dnest_gsl_r;

typedef struct 
{
  double value;
  double tiebreaker;
}LikelihoodType;

typedef struct
{
  LikelihoodType log_likelihood;
  double log_X;
  unsigned long long int visits, exceeds;
  unsigned long long int accepts, tries;
}Level;

extern DNestOptions options;
extern char options_file[STR_MAX_LENGTH];

typedef struct
{
  int id;
  void *addr;
  char tag[50];
  int isset;
}DNestPARDICT;

extern void *particles;
extern int dnest_size_of_modeltype;
extern int particle_offset_size, particle_offset_double;

// sampler
extern bool save_to_disk;
extern unsigned int num_threads;
extern double compression;
extern unsigned int regularisation;

extern void *particles;
extern LikelihoodType *log_likelihoods;
extern unsigned int *level_assignments;

// number account of unaccepted times
extern unsigned int *account_unaccepts;

extern int size_levels, size_levels_combine;  
extern Level *levels;
extern Level *copies_of_levels, *levels_combine;
extern LikelihoodType *all_above;
extern unsigned int count_saves, count_start, num_saves, num_saves_restart;
extern unsigned long long int count_mcmc_steps;
extern LikelihoodType *above;
extern unsigned int size_above, size_all_above, size_all_above_max;

extern int dnest_flag_restart, dnest_flag_postprc, dnest_flag_sample_info, dnest_flag_limits;
extern double dnest_post_temp, dnest_compression;
extern char file_restart[STR_MAX_LENGTH], file_save_restart[STR_MAX_LENGTH];

extern double post_logz;
extern int dnest_num_params;
extern char dnest_sample_postfix[STR_MAX_LENGTH], dnest_sample_tag[STR_MAX_LENGTH], dnest_sample_dir[STR_MAX_LENGTH];
extern double *dnest_param_range, *dnest_prior_info;
extern int *dnest_prior_type;
extern void *dnest_args;

//the limits of parameters for each level;
extern double *limits, *copies_of_limits;

extern int dnest_which_particle_update; // which particle to be updated
extern int dnest_which_level_update;    // which level to be updated;
extern int dnest_thistask, dnest_totaltask;
extern int *dnest_perturb_accept;
extern int dnest_root;
//***********************************************
/*                  functions                  */
extern double mod(double y, double x);
extern void wrap_limit(double *x, double min, double max);
extern int mod_int(int y, int x);
extern int dnest_cmp(const void *pa, const void *pb);

extern void options_load(char *optfile, DNestOptions *opts);
extern void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, 
           double *param_range, int *prior_type, double *prior_info, 
           char *sample_dir, char *optfile, DNestOptions *opts, void *args);
extern void finalise();

extern void dnest_run();
extern void dnest_mcmc_run();
extern void update_particle(unsigned int which);
extern void update_level_assignment(unsigned int which);
extern double log_push(unsigned int which_level);
extern bool enough_levels(Level *l, int size_l);
extern void do_bookkeeping();
extern void save_levels();
extern void save_particle();
extern void save_limits();
extern void kill_lagging_particles();
extern void renormalise_visits();
extern void recalculate_log_X();
extern void dnest_postprocess(double temperature,char *optfile, DNestOptions *opts);
extern void postprocess(double temperature);
extern void initialize_output_file();
extern void close_output_file();
extern void dnest_save_restart();
extern void dnest_restart();
extern void dnest_restart_action(int iflag);
extern void dnest_accept_action();
extern void dnest_kill_action(int i, int i_copy);
extern void dnest_from_prior(void *model);
extern double dnest_perturb(void *model);
extern double dnest_perturb_limit(void *model);
extern void dnest_print_particle(FILE *fp, const void *model);
extern void dnest_read_particle(FILE *fp, void *model);
extern void dnest_check_directory(char *sample_dir);
/*=====================================================*/
// users responsible for following functions
extern void (*print_particle)(FILE *fp, const void *model);
extern void (*read_particle)(FILE *fp, void *model);
extern void (*from_prior)(void *model);
extern double (*log_likelihoods_cal)(const void *model);
extern double (*log_likelihoods_cal_initial)(const void *model);
extern double (*log_likelihoods_cal_restart)(const void *model);
extern double (*perturb)(void *model);
extern void (*restart_action)(int iflag);
extern void (*accept_action)();
extern void (*kill_action)(int i, int i_copy);
/*=====================================================*/

#ifdef __cplusplus
}
#endif

#endif

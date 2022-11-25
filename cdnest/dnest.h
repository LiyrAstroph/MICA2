/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */
#ifndef _DNEST_H
#define _DNEST_H

#ifdef __cplusplus
extern "C" {
#endif

#define DNEST_MAJOR_VERSION 0  // Dec 2, 2018
#define DNEST_MINOR_VERSION 1
#define DNEST_PATCH_VERSION 0

#define STR_MAX_LENGTH (100)
#define BUF_MAX_LENGTH (200)
#define LEVEL_NUM_MAX (1000)

typedef struct
{
  void (*from_prior)(void *model);
  double (*log_likelihoods_cal)(const void *model);
  double (*log_likelihoods_cal_initial)(const void *model);
  double (*log_likelihoods_cal_restart)(const void *model);
  double (*perturb)(void *model);
  void (*print_particle)(FILE *fp, const void *model);
  void (*read_particle)(FILE *fp, void *model);
  void (*restart_action)(int iflag);
  void (*accept_action)();
  void (*kill_action)(int i, int i_copy);
}DNestFptrSet;

// struct for options
typedef struct
{
  unsigned int num_particles;
  unsigned int new_level_interval;
  unsigned int save_interval;
  unsigned int thread_steps;
  unsigned int max_num_levels;
  double lam, beta, max_ptol;
  unsigned int max_num_saves;
  double thread_steps_factor, new_level_interval_factor, save_interval_factor;

  char sample_file[STR_MAX_LENGTH];
  char sample_info_file[STR_MAX_LENGTH];
  char levels_file[STR_MAX_LENGTH];
  char sampler_state_file[STR_MAX_LENGTH];
  char posterior_sample_file[STR_MAX_LENGTH];
  char posterior_sample_info_file[STR_MAX_LENGTH];
  char limits_file[STR_MAX_LENGTH];
}DNestOptions;

extern double dnest_randh();
extern double dnest_rand();
extern double dnest_randn();
extern int dnest_rand_int(int size);
extern void dnest_wrap(double *x, double min, double max);

extern int dnest_get_size_levels();
extern int dnest_get_which_level_update();
extern int dnest_get_which_particle_update();
extern void dnest_get_posterior_sample_file(char *fname);
extern int dnest_check_version(char *verion_str);
extern unsigned int dnest_get_which_num_saves();
extern unsigned int dnest_get_count_saves();
extern unsigned long long int dnest_get_count_mcmc_steps();
extern void dnest_get_limit(int ilevel, int jparam, double *limit1, double *limit2);
extern void dnest_check_fptrset(DNestFptrSet *fptrset);
extern DNestFptrSet * dnest_malloc_fptrset();
extern void dnest_free_fptrset(DNestFptrSet * fptrset);

extern double dnest(int argc, char **argv, DNestFptrSet *fptrset,  int num_params,  
             double *param_range, int *prior_type, double *prior_info, 
             char *sample_dir, char *optfile, DNestOptions *opts, void *args);
             
#ifdef __cplusplus
}
#endif

#endif
/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"
#include "dnestvars.h"

#include "dnest_con.h"
#include "proto.h"

DNestFptrSet *fptrset_con;

double dnest_con(int argc, char **argv)
{
  int i;
  double logz;

  num_params = num_params_var;

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  fptrset_con = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_con->from_prior = from_prior_con;
  fptrset_con->perturb = perturb_con;
  fptrset_con->print_particle = print_particle_con;
  fptrset_con->restart_action = restart_action_con;
  fptrset_con->log_likelihoods_cal = log_likelihoods_cal_con;
  fptrset_con->log_likelihoods_cal_initial = log_likelihoods_cal_initial_con;
  fptrset_con->log_likelihoods_cal_restart = log_likelihoods_cal_restart_con;
  fptrset_con->accept_action = accept_action_con;
  fptrset_con->kill_action = kill_action_con;
  fptrset_con->read_particle = read_particle_con;
  
  set_par_range_con();

  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err == 0)
  {
    for(i=0; i<num_params_var; i+=3)
    {
      par_fix[i] = 1;
      par_fix_val[i] = log(1.0);
    }
  }

  logz = dnest(argc, argv, fptrset_con, num_params, dnest_options_file);

  //free memory
  dnest_free_fptrset(fptrset_con);
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);

  return logz;
}

void set_par_range_con()
{
  int i, j;

  i = 0;
  
  for(j = 0; j<num_params_var; j+=3)
  {
    /* systematic error of continuum */
    par_range_model[i][0] = var_range_model[0][0];
    par_range_model[i++][1] = var_range_model[0][1];

    /* sigma */
    par_range_model[i][0] = var_range_model[1][0];
    par_range_model[i++][1] = var_range_model[1][1];

    /* tau */
    par_range_model[i][0] = var_range_model[2][0];
    par_range_model[i++][1] = var_range_model[2][1];

  }
}

void from_prior_con(void *model)
{
  int i;
  double *pm = (double *)model;
  
  for(i=0; i<num_params; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  return;
}

void print_particle_con(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}


double log_likelihoods_cal_con(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}


double log_likelihoods_cal_initial_con(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}

double log_likelihoods_cal_restart_con(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}

double perturb_con(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width;
  int which, which_level;
  
  /* sample variability parameters more frequently */
  do
  {
    which = dnest_rand_int(num_params);
    
  }while(par_fix[which] == 1);

  /* level-dependent width */
  which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels - 10)?(size_levels-10):which_level_update;
  if( which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }

  pm[which] += dnest_randh() * width;
  dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  
  return logH;
}

void restart_action_con(int iflag)
{
  return;
}

void read_particle_con(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "Error: Cannot read sample file.\n");
      exit(0);
    }
  }

  return;
}

void kill_action_con(int i, int i_copy)
{
  return;
}

void accept_action_con()
{
  return;
}

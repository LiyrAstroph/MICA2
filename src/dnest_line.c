/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "allvars.h"
#include "dnestvars.h"

#include "dnest_line.h"

#include "proto.h"

/* function set for DNest */
DNestFptrSet *fptrset_line;

int dnest_line(int argc, char **argv)
{
  int i;

  /* number of parameters for line */
  if(parset.flag_uniform_tranfuns == 1)
  {
    num_params_line = 4 * nlset_max;   /* each line has a Gussian */
  }
  else
  {
    num_params_line = 0;
    for(i=0; i<nset; i++)
    {
      num_params_line += 4 * dataset[i].nlset; 
    }
  }

  num_params = num_params_var + num_params_line;

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  fptrset_line = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_line->from_prior = from_prior_line;
  fptrset_line->perturb = perturb_line;
  fptrset_line->print_particle = print_particle_line;
  fptrset_line->restart_action = restart_action_line;

  fptrset_line->log_likelihoods_cal = log_likelihoods_cal_line;
  fptrset_line->log_likelihoods_cal_initial = log_likelihoods_cal_initial_line;
  fptrset_line->log_likelihoods_cal_restart = log_likelihoods_cal_restart_line;
  
  set_par_range_line();

  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  /* fix systematic error of continuum */
  if(parset.flag_line_sys_err != 0)
  {
    for(i=0; i<num_params_var; i+=3)
    {
      par_fix[i] = 1;
      par_fix_val[i] = log(1.0);
    }
  }

  /* fix systematic error of line */
  if(parset.flag_line_sys_err != 0)
  {
    for(i=0; i<num_params_line; i+=4)
    {
      par_fix[i+num_params_var] = 1;
      par_fix_val[i+num_params_var] = log(1.0);
    }
  }

  dnest(argc, argv, fptrset_line, num_params, dnest_options_file);

  //free memory
  dnest_free_fptrset(fptrset_line);
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);

  return 0;
}

void set_par_range_line()
{
  int i, j;

  i = 0;
  
  for(j = 0; j<num_params_var; j+=3)
  {
    /* systematic error of continuum */
    par_range_model[i][0] = var_param[j] - 5.0 * var_param_std[j];
    par_range_model[i][1] = var_param[j] + 5.0 * var_param_std[j];

    /* systematic error */
    par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[0][0]);
    par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[0][1]);
    i++;


    /* sigma */
    par_range_model[i][0] = var_param[j+1] - 5.0 * var_param_std[j+1];
    par_range_model[i][1] = var_param[j+1] + 5.0 * var_param_std[j+1];

    par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[1][0]);
    par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[1][1]);
    i++;
   
    /* tau */
    par_range_model[i][0] = var_param[j+2] - 5.0 * var_param_std[j+2];
    par_range_model[i][1] = var_param[j+2] + 5.0 * var_param_std[j+2];

    par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[2][0]);
    par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[2][1]);
    i++;
  }

  for(j = 0; j < num_params_line; j+= 4)
  {
    /* systematic error of line */
    par_range_model[i][0] = line_range_model[0][0];
    par_range_model[i++][1] = line_range_model[0][1];

    /* amplitude of gaussian */
    par_range_model[i][0] = line_range_model[1][0];
    par_range_model[i++][1] = line_range_model[1][1];

    /* center of gaussian */
    par_range_model[i][0] = line_range_model[2][0];
    par_range_model[i++][1] = line_range_model[2][1];

    /* sigma of gaussian */
    par_range_model[i][0] = line_range_model[3][0];
    par_range_model[i++][1] = line_range_model[3][1];
  }
  
  return;
}


void from_prior_line(void *model)
{
  int i;
  double *pm = (double *)model;
  
  for(i=0; i<num_params_var; i++)
  {
    pm[i] = dnest_randn()*var_param_std[i] + var_param[i];
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  
  for(i=num_params_var; i<num_params; i++)
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

void print_particle_line(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}


double log_likelihoods_cal_line(const void *model)
{
  double logL;
  logL = prob_line_variability(model);
  return logL;
}


double log_likelihoods_cal_initial_line(const void *model)
{
  double logL;
  logL = prob_line_variability(model);
  return logL;
}

double log_likelihoods_cal_restart_line(const void *model)
{
  double logL;
  logL = prob_line_variability(model);
  return logL;
}

double perturb_line(void *model)
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

  if(which < num_params_var)
  {
    logH -= (-0.5*pow((pm[which]-var_param[which])/var_param_std[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which]-var_param[which])/var_param_std[which], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  
  return logH;
}


int get_num_params_line()
{
  return num_params;
}

void restart_action_line(int iflag)
{
  return;
}

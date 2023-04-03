/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "allvars.h"
#include "dnest.h"

#include "dnest_pmap.h"

#include "proto.h"

/* function set for DNest */
DNestFptrSet *fptrset_pmap;

double dnest_pmap(int argc, char **argv)
{
  int i, j, k, ic;
  double logz, dlag;

  fptrset_pmap = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_pmap->from_prior = from_prior_pmap;
  fptrset_pmap->perturb = perturb_pmap;
  fptrset_pmap->print_particle = print_particle_pmap;
  fptrset_pmap->restart_action = restart_action_pmap;

  fptrset_pmap->log_likelihoods_cal = log_likelihoods_cal_pmap;
  fptrset_pmap->log_likelihoods_cal_initial = log_likelihoods_cal_initial_pmap;
  fptrset_pmap->log_likelihoods_cal_restart = log_likelihoods_cal_restart_pmap;
  
  fptrset_pmap->accept_action = accept_action_pmap;
  fptrset_pmap->kill_action = kill_action_pmap;
  fptrset_pmap->read_particle = read_particle_pmap;

  /* number of parameters for line */
  if(parset.flag_uniform_tranfuns == 1)
  {
    num_params_line = (1+3*num_gaussian) * nlset_max;   /* each line has the same Gussians */
  }
  else
  {
    num_params_line = 0;
    for(i=0; i<nset; i++)
    {
      num_params_line += (1+3*num_gaussian) * dataset[i].nlset; 
    }
  }

  num_params = num_params_var + num_params_line;

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  model_trans = (void *)malloc(num_params * sizeof(double));

  set_par_range_pmap();
  set_idx_line_pm();

  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err == 0)
  {
    for(i=0; i<num_params_var; i+=3)
    {
      par_fix[i] = 1;
      par_fix_val[i] = log(1.0);
    }
  }

  /* fix systematic error of line */
  if(parset.flag_line_sys_err == 0)
  {
    for(i=0; i<num_params_line; i+=(1+3*num_gaussian))
    {
      par_fix[i+num_params_var] = 1;
      par_fix_val[i+num_params_var] = log(1.0);
    }
  }
  
  print_para_names_pmap();

  logz = dnest(argc, argv, fptrset_pmap, num_params, NULL, NULL, NULL, "data/", dnest_options_file, NULL, NULL);

  //free memory
  dnest_free_fptrset(fptrset_pmap);
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);
  
  free(par_fix);
  free(par_fix_val);

  free(model_trans);

  return logz;
}

void set_par_range_pmap()
{
  int i, j, k;

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

  for(j = 0; j < num_params_line; j+= (1+3*num_gaussian))
  {
    /* systematic error of line */
    par_range_model[i][0] = line_range_model[0][0];
    par_range_model[i++][1] = line_range_model[0][1];

    for(k=0; k<num_gaussian; k++)
    {
      /* amplitude of gaussian */
      if(k==0)
      {
        par_range_model[i][0] = line_range_model[1][0];
        par_range_model[i++][1] = line_range_model[1][1];
      }
      else 
      {
        par_range_model[i][0] = resp_ratio_prior[(k-1)*2+0];
        par_range_model[i++][1] = resp_ratio_prior[(k-1)*2+1];
      }

      /* center of gaussian */
      par_range_model[i][0] = lag_prior[k*2+0]; 
      par_range_model[i++][1] = lag_prior[k*2+1];

      /* sigma of gaussian */
      par_range_model[i][0] = line_range_model[3][0];
      par_range_model[i++][1] = line_range_model[3][1];
    }
  }
  
  return;
}

void print_para_names_pmap()
{
  if(thistask != roottask)
    return;

  int i, j, k;
  char fname[MICA_MAX_STR_LENGTH], fstr[25];

  FILE *fp;

  sprintf(fname, "%s/data/%s_%d", parset.file_dir, "para_names_line.txt", num_gaussian);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  i=-1;
  for(j=0; j<num_params_var; j+=3)
  {
    i++;
    fprintf(fp, "%2d %-25s LOGUNI %10.6f %10.6f %4d %15.6e\n", i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], 
            par_fix[i], par_fix_val[i]);

    i++;
    fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, "sigmad", par_range_model[i][0], par_range_model[i][1], 
            par_fix[i], par_fix_val[i]);

    i++;
    fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, "taud", par_range_model[i][0], par_range_model[i][1], 
            par_fix[i], par_fix_val[i]);
  }

  for(j = 0; j < num_params_line; j+= (1+3*num_gaussian))
  {
    i++;
    fprintf(fp, "%2d %-25s LOGUNI %10.6f %10.6f %4d %15.6e\n", i, "sys_err_line", par_range_model[i][0], par_range_model[i][1], 
            par_fix[i], par_fix_val[i]);

    for(k=0; k<num_gaussian; k++)
    {
      i++;
      if(k == 0)
      {
        sprintf(fstr, "%d-th Gaussian %s", k, "amplitude");
        fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
                par_fix[i], par_fix_val[i]);
      }
      else 
      {
        sprintf(fstr, "%d-th Gaussian %s", k, "resp ratio");
        fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
                par_fix[i], par_fix_val[i]);
      }

      i++;
      sprintf(fstr, "%d-th Gaussian %s", k, "center");
      fprintf(fp, "%2d %-25s UNI    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
              par_fix[i], par_fix_val[i]);

      i++;
      sprintf(fstr, "%d-th Gaussian %s", k, "sigma");
      fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
              par_fix[i], par_fix_val[i]);
    }
  }

  fclose(fp);
  return;
}

void from_prior_pmap(void *model)
{
  int i, j, ic;
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

  /* constrain Gaussian widths */
  if(parset.flag_lag_posivity != 0)
  {
    for(i=0; i<num_params_line; i+= 1+3*num_gaussian)
    {
      ic = num_params_var + i + 2;
      j = 0;
      while(pm[ic-1] - width_factor*pm[ic] < 0.0)
      {
        pm[ic-1] = par_range_model[ic-1][0] + dnest_rand()*(par_range_model[ic-1][1] - par_range_model[ic-1][0]);
        pm[ic] = par_range_model[ic][0] + dnest_rand()*(par_range_model[ic][1] - par_range_model[ic][0]);
        j+1;
        if(j>100)
        {
          pm[ic] = pm[ic-1]/width_factor; 
          break;
        }
      }
    }
  }
  
  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  return;
}

void print_particle_pmap(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}


double log_likelihoods_cal_pmap(const void *model)
{
  double logL;
  logL = prob_line_variability4_pmap(model);
  return logL;
}


double log_likelihoods_cal_initial_pmap(const void *model)
{
  double logL;
  logL = prob_line_variability4_pmap(model);
  return logL;
}

double log_likelihoods_cal_restart_pmap(const void *model)
{
  double logL;
  logL = prob_line_variability4_pmap(model);
  return logL;
}


double perturb_pmap(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width;
  int which, which_level, igau, size_levels;
  
  /* sample variability parameters more frequently */
  do
  {
    which = dnest_rand_int(num_params);
    
  }while(par_fix[which] == 1);

  /* level-dependent width */
  which_level_update = dnest_get_which_level_update();
  size_levels = dnest_get_size_levels();
  
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
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);

    if(parset.flag_lag_posivity)
    {
      int idx;
      idx = check_gauss_positivity(which);
      if(idx == 1 && (pm[which] - width_factor*exp(pm[which+1]) < 0.0))
      {
        logH = -DBL_MAX; /* give a huge penalty */
      }
      else if(idx == 2 && (pm[which-1] - width_factor*exp(pm[which]) < 0.0))
      {
        logH = -DBL_MAX; /* give a huge penalty */
      }
    }
  }
  
  return logH;
}

inline int check_gauss_positivity(int which)
{
  return ((which - num_params_var)%(1+3*num_gaussian)-1)%3;
}

int get_num_params_pmap()
{
  return num_params;
}

void restart_action_pmap(int iflag)
{
  return;
}

void read_particle_pmap(FILE *fp, void *model)
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

void kill_action_pmap(int i, int i_copy)
{
  return;
}

void accept_action_pmap()
{
  return;
}
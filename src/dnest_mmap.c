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
#include <string.h>

#include "allvars.h"
#include "dnest.h"

#include "dnest_mmap.h"

#include "proto.h"

/* function set for DNest */
DNestFptrSet *fptrset_mmap;

double dnest_mmap(int argc, char **argv)
{
  int i, j, k, ic, idx;
  double logz, dlag;
  char dnest_data_dir[MICA_MAX_STR_LENGTH];

  fptrset_mmap = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_mmap->from_prior = from_prior_mmap;
  
  if(type_lag_prior_pr == 0)
    fptrset_mmap->perturb = perturb_mmap_prior0;
  else /* prior 1 or 2 or 3 or 4 */
    fptrset_mmap->perturb = perturb_mmap_prior1;

  fptrset_mmap->print_particle = print_particle_mmap;
  fptrset_mmap->restart_action = restart_action_mmap;

  fptrset_mmap->log_likelihoods_cal = log_likelihoods_cal_mmap;
  fptrset_mmap->log_likelihoods_cal_initial = log_likelihoods_cal_initial_mmap;
  fptrset_mmap->log_likelihoods_cal_restart = log_likelihoods_cal_restart_mmap;
  
  fptrset_mmap->accept_action = accept_action_mmap;
  fptrset_mmap->kill_action = kill_action_mmap;
  fptrset_mmap->read_particle = read_particle_mmap;

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

  set_par_range_line();
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

  /* fix Gaussian lags if type_lag_prior == 2 */
  if(parset.type_lag_prior >= 2 && parset.type_lag_prior <=3)
  {
    dlag = (parset.lag_limit_upper - parset.lag_limit_low)/(num_gaussian-1);
    if(parset.flag_uniform_tranfuns != 1)
    {
      ic = num_params_var;
      for(j=0; j<nset; j++)
      {
        for(k=0; k<dataset[j].nlset; k++)
        {
          ic += 1;  /* error parameter of each line */
          for(i=0; i<num_gaussian; i++)
          {
            par_fix[ic + 1] = 1;  
            par_fix_val[ic + 1] = parset.lag_limit_low + i*dlag;
            ic += 3; 
          }
        }
      }
    }
    else   /* each line has the same Gussians */
    {
      for(j=0; j<nlset_max; j++)
      {
        for(i=0; i<num_gaussian; i++)
        {
          ic = num_params_var + j*(1+3*num_gaussian) + 1 + i*3 + 1;
          par_fix[ic] = 1;
          par_fix_val[ic] = parset.lag_limit_low + i*dlag;
        }
      }
    }
  }

  /* center and width are uniform */
  if(parset.flag_uniform_tranfuns == 2)
  {
    /* start from the second set */
    for(i=1; i<nset; i++)
    {
      for(j=0; j<dataset[i].nlset; j++)
      {
        idx = idx_line_pm[i][j];
        for(k=0; k<num_gaussian; k++)
        {
          par_fix[idx+1 + k*3 + 1] = 1;
          par_fix_val[idx+1 + k*3 + 1] = 99.9; /* center, will be changed to the same as 1st set in perturb */

          par_fix[idx+1 + k*3 + 2] = 1;
          par_fix_val[idx+1 + k*3 + 2] = 99.9; /* width, will be changed to the same as 1st set in perturb */
        }
      }
    }
  }

  if(flag_load_prior == 1)
  {
    load_par_names(prior_file);
    
    /* cope with drw parameters */
    set_drw_par_range_load();

    for(i=0; i<num_params; i++)
    {
      MPI_Bcast(par_range_model[i], 2, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    }
    MPI_Bcast(par_fix, num_params, MPI_INT, roottask, MPI_COMM_WORLD);
    MPI_Bcast(par_fix_val, num_params, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }

  print_para_names_line();
  
  if(flag_para_name != 1)
  {
    strcpy(dnest_data_dir, parset.file_dir);
    strcat(dnest_data_dir, "/data/");
    logz = dnest(argc, argv, fptrset_mmap, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);
  }
  
  //free memory
  dnest_free_fptrset(fptrset_mmap);

  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);

  return logz;
}


void from_prior_mmap(void *model)
{
  int i, j, ic;
  double *pm = (double *)model;
  
  for(i=0; i<num_params_var; i++)
  {
    pm[i] = dnest_randn()*(par_range_model[i][1] - par_range_model[i][0]);
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  
  for(i=num_params_var; i<num_params; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  /* sort component centers */
  if(type_lag_prior_pr == 0)
  {
    double *centers;
    centers = malloc(num_gaussian * sizeof(double));

    for(i=0; i<num_params_line; i+= 1+3*num_gaussian)
    {
      ic = num_params_var + i + 1;
      for(j=0; j<num_gaussian; j++)
      {
        centers[j] = pm[ic+j*3 + 1];
      }
      qsort(centers, num_gaussian, sizeof(double), mica_cmp);

      for(j=0; j<num_gaussian; j++)
      {
        pm[ic + j*3 + 1] = centers[j];
      }
    }
    free(centers);
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

  if(parset.flag_uniform_tranfuns == 2)
  {
    set_center_width(model);
  }

  return;
}

void print_particle_mmap(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}

double log_likelihoods_cal_mmap(const void *model)
{
  double logL;
  logL = prob_line_variability_mmap(model);
  return logL;
}

double log_likelihoods_cal_initial_mmap(const void *model)
{
  double logL;
  logL = prob_line_variability_mmap(model);
  return logL;
}

double log_likelihoods_cal_restart_mmap(const void *model)
{
  double logL;
  logL = prob_line_variability_mmap(model);
  return logL;
}

/* need to check that lag0 < lag1 < lag2... */
double perturb_mmap_prior0(void *model)
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
    if(check_gauss_center(which, &igau) == 1) // Gaussian center is perturbed
    {
      if(igau == 0)
      {
        dnest_wrap(&pm[which], par_range_model[which][0], pm[which+3]);
      }
      else if(igau < num_gaussian - 1)
      {
        dnest_wrap(&pm[which], pm[which-3], pm[which+3]);        
      }
      else
      {
        dnest_wrap(&pm[which], pm[which-3], par_range_model[which][1]);
      }
    }
    else
    {
      dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    } 

    if(parset.flag_lag_posivity)
    {
      logH += check_positivity(model, which);
    }
    if(parset.flag_gap == 1)
    {
      logH += check_gap(model, which);
    }
  }

  if(parset.flag_uniform_tranfuns == 2)
  {
    set_center_width(model);
  }
  
  return logH;
}

/* need not to check lag orders */
double perturb_mmap_prior1(void *model)
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
  
  which_level = which_level_update-5 > (size_levels - 15)?(size_levels-15):which_level_update-5;
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
      logH += check_positivity(model, which);
    }
    if(parset.flag_gap == 1)
    {
      logH += check_gap(model, which);
    }
  }

  if(parset.flag_uniform_tranfuns == 2)
  {
    set_center_width(model);
  }
  
  return logH;
}

int get_num_params_mmap()
{
  return num_params;
}

void restart_action_mmap(int iflag)
{
  return;
}

void read_particle_mmap(FILE *fp, void *model)
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

void kill_action_mmap(int i, int i_copy)
{
  return;
}

void accept_action_mmap()
{
  return;
}
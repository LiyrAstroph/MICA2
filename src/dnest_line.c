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

#include "dnest_line.h"

#include "proto.h"

/* function set for DNest */
DNestFptrSet *fptrset_line;

double dnest_line(int argc, char **argv)
{
  int i, j, k, ic;
  double logz, dlag;
  char *dnest_data_dir[MICA_MAX_STR_LENGTH];

  fptrset_line = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_line->from_prior = from_prior_line;
  if(type_lag_prior_pr == 0)
    fptrset_line->perturb = perturb_line_prior0;
  else /* prior 1 or 2 or 3 or 4 */
    fptrset_line->perturb = perturb_line_prior1;
  fptrset_line->print_particle = print_particle_line;
  fptrset_line->restart_action = restart_action_line;

  fptrset_line->log_likelihoods_cal = log_likelihoods_cal_line;
  fptrset_line->log_likelihoods_cal_initial = log_likelihoods_cal_initial_line;
  fptrset_line->log_likelihoods_cal_restart = log_likelihoods_cal_restart_line;
  
  fptrset_line->accept_action = accept_action_line;
  fptrset_line->kill_action = kill_action_line;
  fptrset_line->read_particle = read_particle_line;

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
    if(parset.flag_uniform_tranfuns == 0)
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

  /* fix tophat widths if type_lag_prior == 3 and type_tf = 1 */
  if(parset.type_tf == 1 && parset.type_lag_prior == 3)
  {
    dlag = (parset.lag_limit_upper - parset.lag_limit_low)/(num_gaussian-1);
    if(parset.flag_uniform_tranfuns == 0)
    {
      ic = num_params_var;
      for(j=0; j<nset; j++)
      {
        for(k=0; k<dataset[j].nlset; k++)
        {
          ic += 1;  /* error parameter of each line */
          for(i=0; i<num_gaussian; i++)
          {
            par_fix[ic + 2] = 1;  
            par_fix_val[ic + 2] = log(dlag/2.0);
            ic += 3; 
          }
        }
      }
    }
    else   /* each line has the same Gussians/tophats */
    {
      for(j=0; j<nlset_max; j++)
      {
        for(i=0; i<num_gaussian; i++)
        {
          ic = num_params_var + j*(1+3*num_gaussian) + 1 + i*3 + 2;
          par_fix[ic] = 1;
          par_fix_val[ic] = log(dlag/2.0);
        }
      }
    }
  }
  
  print_para_names_line();
  
  strcpy(dnest_data_dir, parset.file_dir);
  strcat(dnest_data_dir, "/data/");
  logz = dnest(argc, argv, fptrset_line, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);

  //free memory
  dnest_free_fptrset(fptrset_line);
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);

  return logz;
}

void set_idx_line_pm() 
{
  int i, j, np;

  if(parset.flag_uniform_tranfuns == 1)
  {
    for(i=0; i<nset; i++)
    {
      for(j=0; j<dataset[i].nlset; j++)
      {
        idx_line_pm[i][j] = (1+num_gaussian*3)*j + num_params_var ;
      }
    }
  }
  else
  {
    np = 0;
    for(i=0; i<nset; i++)
    {
      for(j=0; j<dataset[i].nlset; j++)
      {
        idx_line_pm[i][j] = (1+num_gaussian*3)*np + num_params_var;
        np++;
      }
    }
  }

  return;
}

void set_par_range_line()
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
      par_range_model[i][0] = line_range_model[1][0];
      par_range_model[i++][1] = line_range_model[1][1];

      /* center of gaussian */
      if(parset.type_lag_prior == 4)
      {
        par_range_model[i][0] = lag_prior[k*2+0]; 
        par_range_model[i++][1] = lag_prior[k*2+1];
      }
      else if(type_lag_prior_pr == 0)
      {
        par_range_model[i][0] = line_range_model[2][0]; 
        par_range_model[i++][1] = line_range_model[2][1];
      }
      else
      {
        par_range_model[i][0] = line_range_model[2][0] + k*(line_range_model[2][1] - line_range_model[2][0])/(num_gaussian);
        par_range_model[i++][1] =line_range_model[2][0] + (k+1)*(line_range_model[2][1] - line_range_model[2][0])/(num_gaussian);
      }

      /* sigma of gaussian */
      if(strlen(parset.str_width_prior) == 0)
      {  
        if(parset.type_lag_prior == 3 && num_gaussian > 1) /* set the range of sigma to (width/2, width) */
        {
          par_range_model[i][0]   = log((line_range_model[2][1] - line_range_model[2][0])/(num_gaussian-1)/2);
          par_range_model[i++][1] = log((line_range_model[2][1] - line_range_model[2][0])/(num_gaussian-1));
        }
        else 
        {
          par_range_model[i][0]   = line_range_model[3][0];
          par_range_model[i++][1] = line_range_model[3][1];
        }
      }
      else  /* use the input width prior */
      {
        par_range_model[i][0]   = width_prior[k*2+0];
        par_range_model[i++][1] = width_prior[k*2+1];
      }
    }
  }
  
  return;
}


void print_para_names_line()
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
      sprintf(fstr, "%d-th Gaussian %s", k, "amplitude");
      if(parset.flag_negative_resp == 0)
      {
        fprintf(fp, "%2d %-25s LOG    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
                par_fix[i], par_fix_val[i]);
      }
      else 
      {
        fprintf(fp, "%2d %-25s UNI    %10.6f %10.6f %4d %15.6e\n", i, fstr, par_range_model[i][0], par_range_model[i][1], 
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

void from_prior_line(void *model)
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

  /* sort Gaussian centers */
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
  logL = prob_line_variability4(model);
  return logL;
}


double log_likelihoods_cal_initial_line(const void *model)
{
  double logL;
  logL = prob_line_variability4(model);
  return logL;
}

double log_likelihoods_cal_restart_line(const void *model)
{
  double logL;
  logL = prob_line_variability4(model);
  return logL;
}

/* need to check that lag0 < lag1 < lag2... */
double perturb_line_prior0(void *model)
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

/* need not to check lag orders */
double perturb_line_prior1(void *model)
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

/* check where a Gaussian center is updated */
int check_gauss_center(int which, int *igau)
{
  int iwhich, idx;

  iwhich = which - num_params_var;

  idx = (iwhich%(1+3*num_gaussian)-1)%3; //which gaussian parameter is updated; each line has (1+3*num_gaussian) parameters

  *igau = (iwhich%(1+3*num_gaussian)-1)/3; //which gaussian is updated; each gaussian has 3 parameters

  if(idx == 1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
  return 0;
}

int get_num_params_line()
{
  return num_params;
}

void restart_action_line(int iflag)
{
  return;
}

void read_particle_line(FILE *fp, void *model)
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

void kill_action_line(int i, int i_copy)
{
  return;
}

void accept_action_line()
{
  return;
}


int mica_cmp(const void * a, const void * b)
{
  return (*(double*)a) >= (*(double*)b)?1:-1;
}
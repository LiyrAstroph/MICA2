/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

#include "allvars.h"
#include "proto.h"

void init()
{
  int i;

  // order of long-term trend in light curve
  nq = parset.flag_trend+1;

  if(parset.flag_uniform_var_params == 1)
  {
    num_params_var = 3;
  }
  else
  {
    num_params_var = 3*nset;
  }

  /* set parameter range of continuum variation  */
  var_range_model[0][0] = log(1.0); // systematic error in continuum
  var_range_model[0][1] = log(1.0+10.0);

  var_range_model[1][0] = log(1.0e-4); // log(sigma)
  var_range_model[1][1] = log(1.0e1); 

  var_range_model[2][0] = log(1.0); // log(tau)
  var_range_model[2][1] = log(tspan_max*2.0*PI); 


  /* set parameter range of line */
  line_range_model[0][0] = log(1.0); // systematic error in line
  line_range_model[0][1] = log(1.0+10.0);
  
  if(parset.flag_negative_resp == 0)
  {
    line_range_model[1][0] = log(1.0e-5); // amplitude of Gaussian
    line_range_model[1][1] = log(1.0e1);
  }
  else
  {
    line_range_model[1][0] = -1.0e1; // amplitude of Gaussian
    line_range_model[1][1] =  1.0e1;
  }

  line_range_model[2][0] = fmax(parset.lag_limit_low, -tspan_max/2.0); // center of Gaussian
  line_range_model[2][1] = fmin(parset.lag_limit_upper, tspan_max/2.0);

  line_range_model[3][0] = log(tcadence_min/10.0); // sigma of Gaussian
  line_range_model[3][1] = log(tspan_max/3.0); //3*sigma < time span
  if(parset.width_limit_low_isset == 1)
  {
    line_range_model[3][0] = log(parset.width_limit_low);
  }
  if(parset.width_limit_upper_isset == 1)
  {
    /* upper limit cannnot be too large */
    line_range_model[3][1] = fmin(log(parset.width_limit_upper), log(tspan_max/3.0));
  }

  /* cope with str_width_prior */
  width_prior = malloc(parset.num_gaussian_upper * 2 * sizeof(double));
  if(strlen(parset.str_width_prior) > 0)
  {
    char *pstr = parset.str_width_prior;
    int j;

    pstr += 1;
    j = 0;
    for(i=0; i<parset.num_gaussian_upper*2-1; i++)
    {
      sscanf(pstr, "%lf", &width_prior[j]);
      j++;

      pstr = strchr(pstr, ':'); /* values are separated by ":" */
      if(pstr!=NULL)
      {
        pstr++;
      }
      else
      {
        if(thistask == 0)
          printf("No enough width priors.\n");
        exit(0);
      }
    }
    sscanf(pstr, "%lf", &width_prior[j]);
    
    /* convert to log scale */
    for(j=0; j<parset.num_gaussian_upper * 2; j++)
    {
      if(width_prior[j] <= 0)
      {
        if(thistask == 0)
          printf("# Error of %d-th value in str_width_prior!\n", j);
        exit(0);
      }
      else 
      {
        width_prior[j] = log(width_prior[j]);
      }
    }
  }
  else 
  {
    /* set default values */
    for(i=0; i<parset.num_gaussian_upper; i++)
    {
      width_prior[i*2+0] = line_range_model[3][0];
      width_prior[i*2+1] = line_range_model[3][1];
    }
  }

  nscale = 3;
  nrec_max *= nscale;

  if(parset.type_tf == 0)
  {
    if(parset.flag_negative_resp == 0)
    {
      Slc = Slc_gauss;
      Sll = Sll_gauss;
      Sll2 = Sll2_gauss;
      Slc_array = Slc_array_gauss;
      Sll_array = Sll_array_gauss;
      Sll2_array = Sll2_array_gauss;

      Slc_single = Slc_single_gauss;
      Sll_single = Sll_single_gauss;
      Sll2_single = Sll2_single_gauss;
      Sll_single2 = Sll_single2_gauss;
    }
    else 
    {
      Slc = Slc_gauss_linear;
      Sll = Sll_gauss_linear;
      Sll2 = Sll2_gauss_linear;
      Slc_array = Slc_array_gauss_linear;
      Sll_array = Sll_array_gauss_linear;
      Sll2_array = Sll2_array_gauss_linear;

      Slc_single = Slc_single_gauss_linear;
      Sll_single = Sll_single_gauss_linear;
      Sll2_single = Sll2_single_gauss_linear;
      Sll_single2 = Sll_single2_gauss_linear;
    }

    width_factor = 3.0;
  }
  else if(parset.type_tf == 1)
  {
    if(parset.flag_negative_resp == 0)
    {
      Slc = Slc_tophat;
      Sll = Sll_tophat;
      Sll2 = Sll2_tophat;
      Slc_array = Slc_array_tophat;
      Sll_array = Sll_array_tophat;
      Sll2_array = Sll2_array_tophat; 

      Slc_single = Slc_single_tophat;
      Sll_single = Sll_single_tophat;
      Sll2_single = Sll2_single_tophat;
      Sll_single2 = Sll_single2_tophat;
    }
    else
    {
      Slc = Slc_tophat_linear;
      Sll = Sll_tophat_linear;
      Sll2 = Sll2_tophat_linear;
      Slc_array = Slc_array_tophat_linear;
      Sll_array = Sll_array_tophat_linear;
      Sll2_array = Sll2_array_tophat_linear; 

      Slc_single = Slc_single_tophat_linear;
      Sll_single = Sll_single_tophat_linear;
      Sll2_single = Sll2_single_tophat_linear;
      Sll_single2 = Sll_single2_tophat_linear;
    }

    width_factor = 1.0;
  }
  else if(parset.type_tf == 2)
  {
    if(parset.flag_negative_resp == 0)
    {
      Slc = Slc_gamma;
      Sll = Sll_gamma;
      Sll2 = Sll2_gamma;
      Slc_array = Slc_array_gamma;
      Sll_array = Sll_array_gamma;
      Sll2_array = Sll2_array_gamma; 

      Slc_single = Slc_single_gamma;
      Sll_single = Sll_single_gamma;
      Sll2_single = Sll2_single_gamma;
      Sll_single2 = Sll_single2_gamma;
    }
    else
    {
      Slc = Slc_gamma_linear;
      Sll = Sll_gamma_linear;
      Sll2 = Sll2_gamma_linear;
      Slc_array = Slc_array_gamma_linear;
      Sll_array = Sll_array_gamma_linear;
      Sll2_array = Sll2_array_gamma_linear; 

      Slc_single = Slc_single_gamma_linear;
      Sll_single = Sll_single_gamma_linear;
      Sll2_single = Sll2_single_gamma_linear;
      Sll_single2 = Sll_single2_gamma_linear;
    }

    width_factor = 1.0;
  }
  else 
  {
    if(parset.flag_negative_resp == 0)
    {
      Slc = Slc_exp;
      Sll = Sll_exp;
      Sll2 = Sll2_exp;
      Slc_array = Slc_array_exp;
      Sll_array = Sll_array_exp;
      Sll2_array = Sll2_array_exp; 

      Slc_single = Slc_single_exp;
      Sll_single = Sll_single_exp;
      Sll2_single = Sll2_single_exp;
      Sll_single2 = Sll_single2_exp;
    }
    else
    {
      Slc = Slc_exp_linear;
      Sll = Sll_exp_linear;
      Sll2 = Sll2_exp_linear;
      Slc_array = Slc_array_exp_linear;
      Sll_array = Sll_array_exp_linear;
      Sll2_array = Sll2_array_exp_linear; 

      Slc_single = Slc_single_exp_linear;
      Sll_single = Sll_single_exp_linear;
      Sll2_single = Sll2_single_exp_linear;
      Sll_single2 = Sll_single2_exp_linear;
    }

    width_factor = 1.0;
  }

  allocate_memory();
}

void allocate_memory()
{
  int i, ns, n, n_max;

  n_max = nall_max>nrec_max?nall_max:nrec_max;

  PNmat = malloc(nall_max*nall_max*sizeof(double));
  PSmat = malloc(nall_max*nall_max*sizeof(double));
  PCmat = malloc(nall_max*nall_max*sizeof(double));
  IPCmat = malloc(nall_max*nall_max*sizeof(double));

  USmat = malloc(nall_max*nrec_max*sizeof(double));
  USmatT = malloc(nall_max*nrec_max*sizeof(double));

  ASmat = malloc(nrec_max*nrec_max*sizeof(double));

  Tmat1 = malloc(nall_max*nall_max*sizeof(double));
  Tmat2 = malloc(nall_max*nall_max*sizeof(double));
  
  ns = ((1+nlset_max)*nq);
  n = 3*nall_max + (nall_max + 2*nrec_max + 1 + ns) * ns;
  workspace = malloc((n)*sizeof(double));
  workspace_ipiv = malloc(n_max*sizeof(int));
  workspace_inv = malloc(3*nall_max*nall_max*sizeof(double));

  var_param = malloc(num_params_var * sizeof(double));
  var_param_std = malloc(num_params_var * sizeof(double));

  for(i=0; i<num_params_var; i++)
  {
    var_param[i] = var_param_std[i] = 0.0;
  }

  /* index for continuum parameters of each dataset */
  idx_con_pm = malloc(nset*sizeof(int));
  if(parset.flag_uniform_var_params == 1)
  {
    for(i=0; i<nset; i++)
      idx_con_pm[i] = 0;
  }
  else
  {
    for(i=0; i<nset; i++)
      idx_con_pm[i] = 3*i;
  }

  /* evidence */
  logz_arr = malloc(parset.num_gaussian_diff * sizeof(double));

  return;
}

void free_memory()
{
  int i, j;

  /* free dataset */
  for(i=0; i<nset; i++)
  {
    free(dataset[i].con.t);
    free(dataset[i].con.f);
    free(dataset[i].con.fe);

    for(j=0; j<dataset[i].nlset; j++)
    {
      free(dataset[i].line[j].t);
      free(dataset[i].line[j].f);
      free(dataset[i].line[j].fe);
    }
    free(dataset[i].line);
  }
  free(dataset);

  /* free alldata */
  for(i=0; i<nset; i++)
  {
    free(alldata[i].t);
    free(alldata[i].f);
    free(alldata[i].fe);
  }
  free(alldata);

  free(PNmat);
  free(PSmat);
  free(PCmat);
  free(IPCmat);

  free(USmat);
  free(USmatT);
  free(ASmat);
  
  free(Tmat1);
  free(Tmat2);

  free(width_prior);

  free(workspace);
  free(workspace_ipiv);
  free(workspace_inv);

  free(logz_arr);
  return;
}

/* write down OPTIONS file of the continuum run for CDNest */
int write_options_con(char *fname)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    printf("Cannot open file %s\n", fname);
    exit(-1);
  }
  
  fprintf(fp, "# File containing parameters for DNest.\n"
              "# Put comments at the top, or at the end of the line.\n"
              "# Do not change the order of lines.\n"
              "# Lines beginning with '#' are regarded as comments.\n"
              "# DNest options for continuum reconstruction\n\n\n");
  
  fprintf(fp, "MaxNumberSaves            2000\n");
  fprintf(fp, "PTol                      0.1\n");
  fprintf(fp, "NumberParticles           2\n");
  fprintf(fp, "NewLevelIntervalFactor    5\n");
  fprintf(fp, "ThreadStepsFactor         20\n\n\n");
  

  fprintf(fp, "# Full options and their default values (if not specified) are:\n"
              "# MaxNumberSaves           2000 #maximum number of saves\n"
              "# PTol                     0.1   #likelihood tolerance in loge\n"
              "# NumberParticles          2  #number of particles\n"
              "# NewLevelIntervalFactor   5  #new level interval\n"
              "# SaveIntervalFactor       5  #particular saving interval\n"
              "# ThreadStepsFactor        20 #thread steps before communications between cores\n"
              "# MaxNumberLevels          0  #maximum number of levels; unlimited for 0\n"
              "# BacktrackingLength       10.0  #backforward tracking length (lambda)\n"
              "# StrengthEqualPush        100.0 #strength to force equal push (beta)");
  fflush(fp);
  fsync(fileno(fp));

  fclose(fp);

  return 0;
}

/* write down OPTIONS file of the line run for CDNest */
int write_options_line(char *fname)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    printf("Cannot open file %s\n", fname);
    exit(-1);
  }
  
  fprintf(fp, "# File containing parameters for DNest.\n"
              "# Put comments at the top, or at the end of the line.\n"
              "# Do not change the order of lines.\n"
              "# Lines beginning with '#' are regarded as comments.\n"
              "# DNest options for continuum-line time lag analysis\n\n\n");
  
  fprintf(fp, "MaxNumberSaves            %d\n", parset.max_num_saves);
  fprintf(fp, "PTol                      %.2f\n", parset.max_ptol);
  fprintf(fp, "NumberParticles           %d\n",parset.num_particles);
  fprintf(fp, "NewLevelIntervalFactor    %.2f\n", parset.new_level_interval_factor);
  fprintf(fp, "SaveIntervalFactor        %.2f\n", parset.save_interval_factor);
  fprintf(fp, "ThreadStepsFactor         %.2f\n", parset.thread_steps_factor);
  fprintf(fp, "MaxNumberLevels           %d\n", parset.max_num_levels);
  fprintf(fp, "BacktrackingLength        %.2f\n", parset.lam);
  fprintf(fp, "StrengthEqualPush         %.2f\n\n\n", parset.beta);         

  fprintf(fp, "# Full options and their default values (if not specified) are:\n"
              "# MaxNumberSaves           2000 #maximum number of saves\n"
              "# PTol                     0.1   #likelihood tolerance in loge\n"
              "# NumberParticles          1  #number of particles\n"
              "# NewLevelIntervalFactor   1  #new level interval\n"
              "# SaveIntervalFactor       1  #particular saving interval\n"
              "# ThreadStepsFactor        1 #thread steps before communications between cores\n"
              "# MaxNumberLevels          0  #maximum number of levels; unlimited for 0\n"
              "# BacktrackingLength       10.0  #backforward tracking length (lambda)\n"
              "# StrengthEqualPush        100.0 #strength to force equal push (beta)");
  fflush(fp);
  fsync(fileno(fp));

  fclose(fp);

  return 0;
}

void set_argv(int flag_pp, int flag_de, int flag_re)
{
  flag_postprc = flag_pp;
  flag_decomp = flag_de;
  flag_restart = flag_re;
  return;
}
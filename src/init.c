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

#include "allvars.h"
#include "proto.h"

void init()
{
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

  line_range_model[1][0] = log(1.0e-5); // amplitude of Gaussian
  line_range_model[1][1] = log(1.0e2);

  line_range_model[2][0] = fmax(parset.lag_limit_low, -tspan_max/2.0); // center of Gaussian
  line_range_model[2][1] = fmin(parset.lag_limit_upper, tspan_max/2.0);

  line_range_model[3][0] = log(tcadence_min/10.0); // sigma of Gaussian
  line_range_model[3][1] = log(tspan_max/3.0); //3*sigma < time span

  allocate_memory();
}

void allocate_memory()
{
  int i, ns, n;

  PNmat = malloc(nall_max*nall_max*sizeof(double));
  PSmat = malloc(nall_max*nall_max*sizeof(double));
  PCmat = malloc(nall_max*nall_max*sizeof(double));
  IPCmat = malloc(nall_max*nall_max*sizeof(double));

  USmat = malloc(nall_max*nall_max*5*sizeof(double));
  USmatT = malloc(nall_max*nall_max*5*sizeof(double));

  ASmat = malloc(nall_max*nall_max*25*sizeof(double));

  Tmat1 = malloc(nall_max*nall_max*sizeof(double));
  Tmat2 = malloc(nall_max*nall_max*sizeof(double));
  
  ns = ((1+nlset_max)*nq);
  n = 10*nall_max + (6*nall_max + 1 + ns) * ns;
  workspace = malloc((n)*sizeof(double));
  workspace_ipiv = malloc(5*nall_max*sizeof(int));
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
              "# MaxNumberSaves           10000 #maximum number of saves\n"
              "# PTol                     0.1   #likelihood tolerance in loge\n"
              "# NumberParticles          1  #number of particles\n"
              "# NewLevelIntervalFactor   2  #new level interval\n"
              "# SaveIntervalFactor       2  #particular saving interval\n"
              "# ThreadStepsFactor        10 #thread steps before communications between cores\n"
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
              "# MaxNumberSaves           10000 #maximum number of saves\n"
              "# PTol                     0.1   #likelihood tolerance in loge\n"
              "# NumberParticles          1  #number of particles\n"
              "# NewLevelIntervalFactor   2  #new level interval\n"
              "# SaveIntervalFactor       2  #particular saving interval\n"
              "# ThreadStepsFactor        10 #thread steps before communications between cores\n"
              "# MaxNumberLevels          0  #maximum number of levels; unlimited for 0\n"
              "# BacktrackingLength       10.0  #backforward tracking length (lambda)\n"
              "# StrengthEqualPush        100.0 #strength to force equal push (beta)");
  fflush(fp);
  fsync(fileno(fp));

  fclose(fp);

  return 0;
}
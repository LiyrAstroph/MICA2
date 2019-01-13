/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void init()
{
  // order of long-term trend in light curve
  nq = 1;

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

  line_range_model[3][0] = log(tcadence_line_min*0.3); // sigma of Gaussian
  line_range_model[3][1] = log(tspan_max);

  allocate_memory();
}

void allocate_memory()
{
  int i;

  PNmat = malloc(nall_max*nall_max*sizeof(double));
  PSmat = malloc(nall_max*nall_max*sizeof(double));
  PCmat = malloc(nall_max*nall_max*sizeof(double));
  IPCmat = malloc(nall_max*nall_max*sizeof(double));

  USmat = malloc(nall_max*nall_max*5*sizeof(double));
  USmatT = malloc(nall_max*nall_max*5*sizeof(double));

  ASmat = malloc(nall_max*nall_max*25*sizeof(double));

  Tmat1 = malloc(nall_max*nall_max*sizeof(double));
  Tmat2 = malloc(nall_max*nall_max*sizeof(double));

  workspace = malloc((50*nall_max)*sizeof(double));
  workspace_ipiv = malloc(nall_max*sizeof(double));

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

  free(logz_arr);
  return;
}
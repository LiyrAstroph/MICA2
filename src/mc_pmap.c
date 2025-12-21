/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_erf.h>

#include "dnest_pmap.h"
#include "allvars.h"
#include "dnest.h"

#include "proto.h"

void mc_pmap(double *logz)
{
  int i, j, argc=0, jzmax=0;
  double logz_max;
  char **argv;

  argv = malloc(11*sizeof(char *));
  for(i=0; i<11; i++)
  {
    argv[i] = malloc(MICA_MAX_STR_LENGTH*sizeof(char));
  }
  /* setup argc and argv */
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restart_dnest1d.txt");

  /* level-dependent sampling */
  {
    strcpy(argv[argc++], "-l");
  }

  if(flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart_dnest1d.txt");
  }

  if(flag_postprc == 1)
    strcpy(argv[argc++], "-p");
  
  // sample tag 
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "1d");

  mc_pmap_init();

  logz_max = -DBL_MAX;
  for(j=0; j<parset.num_gaussian_diff; j++)
  {
    num_gaussian = parset.num_gaussian_low + j;
    
    if(thistask == roottask)
    {
      printf("# number of components: %d\n", num_gaussian);
    }
    
    sprintf(postfix, "_%d", num_gaussian);
    strcpy(argv[argc], "-x");
    strcpy(argv[argc+1], postfix);

    logz[j] = dnest_pmap(argc+2, argv);
    if(logz_max < logz[j])
    {
      logz_max = logz[j];
      jzmax = j;
    }
    
    if(flag_para_name != 1 && flag_postsample != 1)
    {
      postprocess_line();

      output_reconstruction_parallel();
      
      if(num_gaussian > 1 && flag_decomp == 1)
      {
        output_decompose_line_parallel();
      }
    }
  }
  
  if(thistask == roottask && flag_para_name != 1)
  {
    char fname[200];
    printf("*****************************************************\n");
    for(j=0; j<parset.num_gaussian_diff; j++)
    {
      printf("number of components: %d, evidence: %f\n", parset.num_gaussian_low + j, logz[j]);
    }

    printf("best number components: %d.\n", parset.num_gaussian_low + jzmax);
    printf("*****************************************************\n");

    FILE *fp=NULL;
    sprintf(fname, "%s/%s", parset.file_dir, "data/evidence.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# number_of_components     evidence log(z)\n");
    for(j=0; j<parset.num_gaussian_diff; j++)
    {
      fprintf(fp, "%d       %f\n", parset.num_gaussian_low + j, logz[j]);
    }
    fclose(fp);
  } 

  mc_pmap_end();

  /* clear up argv */
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

int mc_pmap_init()
{
  int i;

  mc_line_init();
  
  resp_ratio_prior = malloc((parset.num_gaussian_upper -1)*2*sizeof(double));
  /* default values */
  if(strlen(parset.str_ratio_prior) == 0)
  {
    for(i=0; i<parset.num_gaussian_upper-1; i++)
    {
      resp_ratio_prior[i*2+0] = log(1.0e-3);
      resp_ratio_prior[i*2+1] = log(1.0);
    }
  }
  else 
  {
    char *pstr = parset.str_ratio_prior;
    int j, jc;

    /* first count the number of ":" */
    pstr += 1;
    jc = 0;
    while(1)
    {
      pstr = strchr(pstr, ':');
      if(pstr != NULL)
      {
        pstr++;
        jc++;
      }
      else 
      {
        break;
      }
    }
    if(jc<1 || jc%2 == 0)
    {
      if(thistask == 0)
        printf("Incorrect StrRatioPrior.\n");
      exit(0);
    }
    
    /* now begin to read off */
    pstr = parset.str_ratio_prior;
    pstr += 1;
    j = 0;
    for(i=0; i<jc*2-1; i++)
    {
      sscanf(pstr, "%lf", &resp_ratio_prior[j]);
      resp_ratio_prior[j] = log(resp_ratio_prior[j]);
      j++;

      pstr = strchr(pstr, ':'); /* values are separated by ":" */
      if(pstr!=NULL)
      {
        pstr++;
      }
    }
    sscanf(pstr, "%lf", &resp_ratio_prior[j]);
    resp_ratio_prior[j] = log(resp_ratio_prior[j]);
    
    /* pad the rest ratios */
    for(i=jc; i<parset.num_gaussian_upper-1; i++)
    {
      resp_ratio_prior[i*2+0] = resp_ratio_prior[(jc-1)*2+0];
      resp_ratio_prior[i*2+1] = resp_ratio_prior[(jc-1)*2+1];
    }
  }

  /* cope with line width, set a lower limit to width of line TF */
  /* the first gaussian is disk RM, width might be small  */
  width_prior[0*2+0] = fmin(width_prior[0*2+0], log(0.01));
  for(i=1; i<parset.num_gaussian_upper; i++)
  {
    width_prior[i*2+0] = fmax(width_prior[i*2+0], log(1.0));
  }

  return 0;
}

int mc_pmap_end()
{
  mc_line_end();
  free(resp_ratio_prior);

  return 0;
}

/*!
 * cope with response ratios 
 * 
 * note that model_in and model_out can be same
 */
void transform_response_ratio(const void *model_in, void *model_out)
{
  int j, k, idx;
  double *pm_in = (double *)model_in;
  double *pm_out = (double *)model_out;
  
  memcpy(model_out, model_in, sizeof(double) * num_params);
  
  idx = num_params_var;
  /* each line has (1+3*num_gaussian) parameters */
  for(j = 0; j < num_params_line; j+= (1+3*num_gaussian))
  {
    for(k=1; k<num_gaussian; k++)
    {
      /* note in log */
      pm_out[idx + j + 1 + k*3] += pm_in[idx + j + 1 + 0*3];
    }
  }

  return;
}
void transform_response_ratio_inplace(const void *model_in)
{
  int j, k, idx;
  double *pm_in = (double *)model_in;
    
  idx = num_params_var;
  /* each line has (1+3*num_gaussian) parameters */
  for(j = 0; j < num_params_line; j+= (1+3*num_gaussian))
  {
    for(k=1; k<num_gaussian; k++)
    {
      /* note in log */
      pm_in[idx + j + 1 + k*3] += pm_in[idx + j + 1 + 0*3];
    }
  }

  return;
}

/*!
 * using semiseparable matrix and partitioning inverse
 * for pmap
 * 
 */
double prob_line_variability4_pmap(const void *model)
{
  double prob = 0.0, prob1, sigma, tau, syserr;
  int i, j, k, m, np, info, sign, *ipiv;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq, *W, *D, *phi, *fe;
  double *fall;
  int nall, nqall, idx;
  double *pm = (double *)model;
  double *work;
  int *narr, nd;

  /* first cope with response ratios */
  transform_response_ratio(model, model_trans);

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max * ((1+nlset_max)*nq);
  yq = y + nall_max;
  Cq = yq + (1+nlset_max)*nq;

  W = Cq + ((1+nlset_max)*nq) * ((1+nlset_max)*nq);
  D = W + nall_max;
  phi = D + nall_max;
  
  fe = phi + nall_max;

  narr = (int *)(fe + nall_max);

  ipiv = workspace_ipiv;
  work = workspace_inv;

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    idx = idx_con_pm[k];
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]) * sqrt(tau);
    syserr = (exp(pm[idx+0])-1.0)*dataset[k].con.error_mean/sigma;

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (1+dataset[k].nlset);

    for(i=0;i<dataset[k].con.n;i++)
    {
      Larr[i*nqall] = 1.0; 
      for(j=1; j<nq; j++)
        Larr[i*nqall+j] = pow(dataset[k].con.t[i], j);

      for(j=nq; j<nqall; j++)
        Larr[i*nqall + j] =  0.0;
    }
    np = dataset[k].con.n;
    for(j=0; j<dataset[k].nlset; j++)
    {
      for(m=0; m<dataset[k].line[j].n; m++)
      {
        for(i=0; i<nqall; i++)
          Larr[(np+m)*nqall + i ]  = 0.0;
        
        Larr[(np+m)*nqall + nq + j*nq + 0] = 1.0;
        for(i=1; i<nq; i++)
          Larr[(np+m)*nqall + nq + j*nq + i] = pow(dataset[k].line[j].t[m], i);
      }
      np += dataset[k].line[j].n;
    }
    
    nd = dataset[k].nlset;
    narr[0] = dataset[k].con.n;
    for(i=1; i<=nd; i++)
      narr[i] = dataset[k].line[i-1].n;

    set_covar_Pmat_data_line_array(model_trans, k);

    for(i=0; i<dataset[k].con.n; i++)
      fe[i] = dataset[k].con.fe[i]/sigma;
      
    inverse_semiseparable_iter(dataset[k].con.t, dataset[k].con.n, 1.0, 1.0/tau, 
                          fe, syserr, W, D, phi, PCmat, nall, narr, nd,
                          &lndet, work, ipiv);

    /* C^-1 and lndet(C) */
    //inverse_symat_partition_iter(PCmat, nall, narr, nd, &lndet, work, ipiv); 
    lndet += 2.0*nall*log(sigma);

    /* calculate L^T*C^-1*L */
    multiply_mat_MN(PCmat, Larr, ybuf, nall, nqall, nall);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall);

    /* calculate L^T*C^-1*y */
    // multiply_matvec(PCmat, fall, nall, ybuf);
    // multiply_mat_MN_transposeA(Larr, ybuf, yq, nqall, 1, nall);
    
    /* L^T*C^-1*y = (C^-1xL)^T x y = ybuf^T x y */
    multiply_matvec_MN_transposeA(ybuf, nall, nqall, fall, yq);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_symat_lndet_sign(Cq, nqall, &lndet_ICq, &info, &sign, ipiv);
    if(info!=0 || sign==-1 )
    {
      prob = -DBL_MAX;
      if(flag_warning)printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
      return prob;
    }
    lndet_ICq += - 2.0*nqall*log(sigma);
    multiply_mat_MN(Cq, yq, ybuf, nqall, 1, nqall);
    
    /* y = y - L x q = -1 * (L x q) + 1 * y */
    // multiply_matvec_MN(Larr, nall, nqall, ybuf, y);
    // for(i=0; i<nall; i++)
    // {
    //   y[i] = fall[i] - y[i];
    // }
    memcpy(y, fall, nall*sizeof(double));
    multiply_matvec_MN_alpha_beta(Larr, nall, nqall, ybuf, y, -1.0, 1.0);

    /* y^T x C^-1 x y */
    multiply_matvec(PCmat, y, nall, ybuf);
    prob1 = -0.5 * cblas_ddot(nall, y, 1, ybuf, 1)/(sigma*sigma);
    
    if(prob1 > 0.0 )  // check if prob is positive
    { 
      prob = -DBL_MAX;
      if(flag_warning)printf("prob >0!\n");
      return prob;
    }
    
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;
  }
  return prob;
}
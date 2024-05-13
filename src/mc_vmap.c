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

#include "dnest_vmap.h"
#include "allvars.h"

#include "proto.h"

void mc_vmap()
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

  mc_vmap_init();

  logz_max = -DBL_MAX;
  for(j=0; j<parset.num_gaussian_diff; j++)
  {
    num_gaussian = parset.num_gaussian_low + j;
    
    if(thistask == roottask)
    {
      printf("# number of Gaussian/Tophat: %d\n", num_gaussian);
    }
    
    sprintf(postfix, "_%d", num_gaussian);
    strcpy(argv[argc], "-x");
    strcpy(argv[argc+1], postfix);

    logz_arr[j] = dnest_vmap(argc+2, argv);
    if(logz_max < logz_arr[j])
    {
      logz_max = logz_arr[j];
      jzmax = j;
    }

    postprocess_vmap();

    output_reconstruction_vmap_parallel();
  }
  
  if(thistask == roottask)
  {
    char fname[200];
    printf("*****************************************************\n");
    for(j=0; j<parset.num_gaussian_diff; j++)
    {
      printf("number of Gaussian/Tophat: %d, evidence: %f\n", parset.num_gaussian_low + j, logz_arr[j]);
    }

    printf("best number of Gaussian/Tophat: %d.\n", parset.num_gaussian_low + jzmax);
    printf("*****************************************************\n");

    FILE *fp;
    sprintf(fname, "%s/%s", parset.file_dir, "data/evidence.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# number_of_gaussians/tophats     evidence\n");
    for(j=0; j<parset.num_gaussian_diff; j++)
    {
      fprintf(fp, "%d       %f\n", parset.num_gaussian_low + j, logz_arr[j]);
    }
    fclose(fp);
  } 

  mc_vmap_end();

  /* clear up argv */
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;

}


/*!
 *  this function does postprocess. 
 */
void postprocess_vmap()
{
  char posterior_sample_file[MICA_MAX_STR_LENGTH];
  double *pm, *pmstd;
  int num_ps, i, j;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);
    
  if(thistask == roottask)
  {
    void *best_model_line;   /*!< best model */
    void *best_model_std_line;  /*!< standard deviation of the best model */
    FILE *fp;

    best_model_line = malloc(size_of_modeltype);
    best_model_std_line = malloc(size_of_modeltype);

    /* get file name of posterior sample file */
    dnest_get_posterior_sample_file(posterior_sample_file);

    /* open file for posterior sample */
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }

    /* read number of points in posterior sample */
    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);
    }

    fclose(fp);

    /* calcaulte mean and standard deviation of posterior samples. */
    pm = (double *)best_model_line;
    pmstd = (double *)best_model_std_line;
    for(j=0; j<num_params; j++)
    {
      pm[j] = pmstd[j] = 0.0;
    }
    for(i=0; i<num_ps; i++)
    {
      for(j =0; j<num_params; j++)
        pm[j] += *((double *)posterior_sample + i*num_params + j );
    }

    for(j=0; j<num_params; j++)
      pm[j] /= num_ps;

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
        pmstd[j] += pow( *((double *)posterior_sample + i*num_params + j ) - pm[j], 2.0 );
    }

    for(j=0; j<num_params; j++)
    {
      if(num_ps > 1)
        pmstd[j] = sqrt(pmstd[j]/(num_ps-1.0));
      else
        pmstd[j] = 0.0;
    }  

    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_line + j), 
                                             *((double *)best_model_std_line+j) ); 

    free(post_model);
    free(posterior_sample);

    free(best_model_line);
    free(best_model_std_line);
  }
  return;
}

int mc_vmap_init()
{
  int i;

  if(thistask == roottask)
  {
    if(strlen(dnest_options_file_line) == 0)
    {
      sprintf(dnest_options_file, "%s/param/OPTIONS1D", parset.file_dir);
      write_options_line(dnest_options_file);
    }
    else 
    {
      strcpy(dnest_options_file, dnest_options_file_line);
      get_num_particles(dnest_options_file);
    }
  }
  MPI_Bcast(dnest_options_file, MICA_MAX_STR_LENGTH, MPI_CHAR, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  idx_line_pm = malloc(nset * sizeof(double *));
  for(i=0; i<nset; i++)
  {
    idx_line_pm[i] = malloc(nlset_max * sizeof(int));
  }

  Smat_lc = malloc(ncon_max * nline_max * sizeof(double));
  Smat_ll = malloc(nline_max * nline_max * sizeof(double));
  
  /* get lag prior from parset.str_lag_prior */
  if(parset.type_lag_prior == 4)
  {
    lag_prior = malloc(parset.num_gaussian_upper * 2 * sizeof(double));
    char *pstr = parset.str_lag_prior;
    int j;

    pstr += 1;
    j = 0;
    for(i=0; i<parset.num_gaussian_upper*2-1; i++)
    {
      sscanf(pstr, "%lf", &lag_prior[j]);
      j++;

      pstr = strchr(pstr, ':'); /* values are separated by ":" */
      if(pstr!=NULL)
      {
        pstr++;
      }
      else
      {
        if(thistask == 0)
          printf("No enough lag priors.");
        exit(0);
      }
    }
    sscanf(pstr, "%lf", &lag_prior[j]);
  }

  return 0;
}

int mc_vmap_end()
{
  int i;

  for(i=0; i<nset; i++)
  {
    free(idx_line_pm[i]);
  }
  free(idx_line_pm);

  free(Smat_lc);
  free(Smat_ll);
  
  if(parset.type_lag_prior == 4)
  {
    free(lag_prior);
  }
  return 0;
}

/*!
 * calulate inverse and lndet simultaneously using direct matrix inverse
 */
double prob_line_variability3_vmap(const void *model)
{
  double prob = 0.0, prob1, sigma, tau;
  int i, j, k, m, np, info, sign, *ipiv;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq;
  double *fall;
  int nall, nqall, idx;
  double *pm = (double *)model;

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max* ((1+nlset_max)*nq);
  yq = y + nall_max;
  Cq = yq + (1+nlset_max)*nq;

  ipiv = workspace_ipiv;

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    idx = idx_con_pm[k];
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]); //* sqrt(tau);

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (dataset[k].nlset);

    /* no continuum data points */
    np = 0;
    for(j=0; j<dataset[k].nlset; j++)
    {
      for(m=0; m<dataset[k].line[j].n; m++)
      {
        for(i=0; i<nqall; i++)
          Larr[(np+m)*nqall + i ]  = 0.0;
        
        Larr[(np+m)*nqall + j*nq + 0] = 1.0;
        for(i=1; i<nq; i++)
          Larr[(np+m)*nqall + j*nq + i] = pow(dataset[k].line[j].t[m], i);
      }
      np += dataset[k].line[j].n;
    }
    
    set_covar_Pmat_data_line_array_vmap(model, k);

    inverse_symat_lndet_sign(PCmat, nall, &lndet, &info, &sign, ipiv); /* calculate C^-1 */
    if(info!=0|| sign==-1)
    {
      prob = -DBL_MAX;
      printf("lndet_C %f %d!\n", lndet, sign);
      return prob;
    }
    lndet += 2.0*nall*log(sigma);

    /* calculate L^T*C^-1*L */
    multiply_mat_MN(PCmat, Larr, ybuf, nall, nqall, nall);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall);

    /* calculate L^T*C^-1*y */
    multiply_matvec(PCmat, fall, nall, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nqall, 1, nall);
    
    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_symat_lndet_sign(Cq, nqall, &lndet_ICq, &info, &sign, ipiv);
    if(info!=0 || sign==-1 )
    {
      prob = -DBL_MAX;
      printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
      return prob;
    }
    lndet_ICq += - 2.0*nqall*log(sigma);
    multiply_mat_MN(Cq, yq, ybuf, nqall, 1, nqall);
  
    multiply_matvec_MN(Larr, nall, nqall, ybuf, y);
    for(i=0; i<nall; i++)
    {
      y[i] = fall[i] - y[i];
    }

    /* y^T x C^-1 x y */
    multiply_matvec(PCmat, y, nall, ybuf);
    prob1 = -0.5 * cblas_ddot(nall, y, 1, ybuf, 1)/(sigma*sigma);
    
    if(prob1 > 0.0 )  // check if prob is positive
    { 
      prob = -DBL_MAX;
      printf("prob >0!\n");
      return prob;
    }
    
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;
  }
  return prob;
}

/*!
 * calulate inverse and lndet simultaneously using matrix partition
 */
double prob_line_variability4_vmap(const void *model)
{
  double prob = 0.0, prob1, sigma, tau;
  int i, j, k, m, np, info, sign, *ipiv;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq;
  double *fall;
  int nall, nqall, idx;
  double *pm = (double *)model;
  double *work;
  int *narr, nd;

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max* ((1+nlset_max)*nq);
  yq = y + nall_max;
  Cq = yq + (1+nlset_max)*nq;

  narr = (int *)(Cq + (1+nlset_max)*nq*(1+nlset_max)*nq);

  ipiv = workspace_ipiv;
  work = workspace_inv;

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    idx = idx_con_pm[k];
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]); //* sqrt(tau);

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (dataset[k].nlset);

    /* no continuum data points */
    np = 0;
    for(j=0; j<dataset[k].nlset; j++)
    {
      for(m=0; m<dataset[k].line[j].n; m++)
      {
        for(i=0; i<nqall; i++)
          Larr[(np+m)*nqall + i ]  = 0.0;
        
        Larr[(np+m)*nqall + j*nq + 0] = 1.0;
        for(i=1; i<nq; i++)
          Larr[(np+m)*nqall + j*nq + i] = pow(dataset[k].line[j].t[m], i);
      }
      np += dataset[k].line[j].n;
    }
    
    nd = dataset[k].nlset - 1;
    narr[0] = dataset[k].line[0].n;
    for(i=1; i<=nd; i++)
      narr[i] = dataset[k].line[i].n;

    set_covar_Pmat_data_line_array_vmap(model, k);

    /* C^-1 */
    inverse_symat_partition_iter(PCmat, nall, narr, nd, &lndet, work, ipiv); 
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
      printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
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
      printf("prob >0!\n");
      return prob;
    }
    
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;
  }
  return prob;
}

/*!
 * this function sets the covariance matrix at data time points for an array of times 
 * for vmap
 * k is the index of dataset
 */
void set_covar_Pmat_data_line_array_vmap(const void *model, int k)
{
  double t1, t2;
  double sigma, tau, syserr_line, error;
  int i, j, m, l, nall, np, npline, idx, nline, nline2;
  double *pm = (double *)model, *tcon, *tline, *tline2;
  
  /* total number of point of this dataset */
  nall = alldata[k].n;

  /* set variability parameters */
  idx = idx_con_pm[k];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]); // * sqrt(tau);

  /* no continuum data point */
  /* line - line */
  np = 0;
  for(j=0; j<dataset[k].nlset; j++)
  {
    idx = idx_line_pm[k][j];
    syserr_line = (exp(pm[idx]) - 1.0) * dataset[k].line[j].error_mean;
    
    tline = dataset[k].line[j].t;
    nline = dataset[k].line[j].n;
    Sll_array(tline, nline, model, k, j, Smat_ll);
    /* line self */
    for(i=0; i<nline; i++)
    {
      for(m=0; m<i; m++)
      {
        PCmat[(np + i)*nall+ (np + m) ] = PCmat[(np+m)*nall+(np+i)] =
        PSmat[(np + i)*nall+ (np + m) ] = PSmat[(np+m)*nall+(np+i)] = Smat_ll[i*nline + m];
      }
      PSmat[(np + i)*nall+(np +i)] = Smat_ll[i*nline+i];
      error = dataset[k].line[j].fe[i] * dataset[k].line[j].fe[i] + syserr_line*syserr_line;
      PCmat[(np + i)*nall+(np +i)] = PSmat[(np + i)*nall+(np +i)] + error/sigma/sigma;
    }

    /* between lines */
    npline = np + nline;
    for(l=j+1; l<dataset[k].nlset; l++)
    {
      tline2 = dataset[k].line[l].t;
      nline2 = dataset[k].line[l].n;
      Sll2_array(tline, nline, tline2, nline2, model, k, j, l, Smat_ll);

      for(i=0; i<nline; i++)
      {
        for(m=0; m<nline2; m++)
        {
          PCmat[(np+i)*nall + npline + m ] = PCmat[ (npline + m) * nall + np + i ] =
          PSmat[(np+i)*nall + npline + m ] = PSmat[ (npline + m) * nall + np + i ] = Smat_ll[i*nline2 + m];
        }
      }
      npline += nline2;
    }
    np += nline;
  }
  
  return;
}

/* 
 * reconstruct line lcs acording to the time span of the continuum in parallel.
 */

void output_reconstruction_vmap_parallel()
{
  int i, j, k, m;
  int num_ps, size_of_modeltype;
  int num_ps_task, num_ps_task_average;
  void *posterior_sample_task, *posterior_sample_all;
  double *post_model, *ps;

  double **tall, **fall, **feall, **feall_max, **fall_best, **fall_std, *yq, **yq_best, **yq_std;
  double **fall_best_buf, **fall_std_buf, **feall_max_buf, **yq_best_buf, **yq_std_buf;
  int **nall, *ntall, np;
  double tspan, tbeg, tend;

  char fname[200];
  FILE *fp, *fp_sample, *fpq;
  
  size_of_modeltype = num_params * sizeof(double);

  if(thistask == roottask)
  {
    printf("start reconstruction...\n");
    
    char posterior_sample_file[MICA_MAX_STR_LENGTH];

    /* get file name of posterior sample file */
    dnest_get_posterior_sample_file(posterior_sample_file);
    
    /* open file for posterior sample */
    fp_sample = fopen(posterior_sample_file, "r");
    if(fp_sample == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    /* read number of points in posterior sample */
    if(fscanf(fp_sample, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }

    posterior_sample_all = (void *)malloc(num_ps * size_of_modeltype);
  
    /* start to read in the posterior sample */
    post_model = (double *)posterior_sample_all;
    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp_sample, "%lf", post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp_sample, "\n");

      post_model += num_params;
    }
    fclose(fp_sample);

    sprintf(fname, "%s/%s%s", parset.file_dir, "data/pall.txt", postfix);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    fprintf(fp, "# %d\n", nset);

    sprintf(fname, "%s/%s%s", parset.file_dir, "data/trend.txt", postfix);
    fpq = fopen(fname, "w");
    if(fpq == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }
  
  MPI_Bcast(&num_ps, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  
  num_ps_task_average  = num_ps/totaltask;
  if(thistask == roottask)
    num_ps_task = num_ps - num_ps_task_average * (totaltask -1);
  else
    num_ps_task = num_ps_task_average;
  
  posterior_sample_task = (void *)malloc(num_ps_task * size_of_modeltype);
  
  if(thistask == roottask)
  {
    /* offset the point of roottask */
    post_model = (double *)posterior_sample_all;
    post_model += (num_ps_task - num_ps_task_average)*num_params;
  }
  /* scatter the posterior sample */
  MPI_Scatter(post_model, num_ps_task_average*num_params, MPI_DOUBLE, 
              posterior_sample_task, num_ps_task_average*num_params, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  
  /* cope with the roottask */
  if(thistask == roottask)
  {
    ps = (double *)posterior_sample_task;
    ps += num_ps_task_average * num_params;
    post_model = (double *)posterior_sample_all;
    memcpy((void *)ps, (void *)post_model, (num_ps_task - num_ps_task_average) * size_of_modeltype);
  }
  
  /* now start to reconstruction */
  nall = malloc(nset*sizeof(int *));
  tall = malloc(nset*sizeof(double *));
  fall = malloc(nset*sizeof(double *));
  feall = malloc(nset*sizeof(double *));
  feall_max = malloc(nset*sizeof(double *));
  fall_best = malloc(nset*sizeof(double *));
  fall_std = malloc(nset*sizeof(double *));
  ntall = malloc(nset*sizeof(int));
  yq_best = malloc(nset*sizeof(double *));
  yq_std = malloc(nset*sizeof(double *));
  
  yq = malloc(nq*(1+nlset_max)*sizeof(double));

  /* first set the number of reconstruction points */
  for(i=0; i<nset; i++)
  {
    nall[i] = malloc((1+nlset_max) * sizeof(int));
    
    nall[i][0] = 0;
    for(j=0; j<dataset[i].nlset; j++)
    {
      /* number < 500 */
      nall[i][0] = fmax(nall[i][0], fmin(dataset[i].line[j].n*nscale, 500));
    }

    for(j=0; j<dataset[i].nlset; j++)
      nall[i][1+j] = nall[i][0];
    
    /* compute total number of points */
    ntall[i] = nall[i][0];
    for(j=0; j<dataset[i].nlset; j++)
      ntall[i] += nall[i][1+j];
  }

  for(i=0; i<nset; i++)
  { 
    tall[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall[i] = (double *)malloc(ntall[i] * sizeof(double));
    feall[i] = (double *)malloc(ntall[i] * sizeof(double));
    feall_max[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_best[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_std[i] = (double *)malloc(ntall[i] * sizeof(double));

    yq_best[i] = (double *)malloc(nq*(1+nlset_max)*sizeof(double));
    yq_std[i] = (double *)malloc(nq*(1+nlset_max)*sizeof(double));
  }

  /* bufs for reduce */
  feall_max_buf = (double **)malloc(nset*sizeof(double *));
  fall_best_buf = (double **)malloc(nset*sizeof(double *));
  fall_std_buf = (double **)malloc(nset*sizeof(double *));
  yq_best_buf = (double **)malloc(nset*sizeof(double *));
  yq_std_buf = (double **)malloc(nset*sizeof(double *));

  for(i=0; i<nset; i++)
  {
    feall_max_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_best_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_std_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
    yq_best_buf[i] = (double *)malloc(nq*(1+nlset_max)*sizeof(double));
    yq_std_buf[i] = (double *)malloc(nq*(1+nlset_max)*sizeof(double));
  }
  
  /* initialize */
  for(i=0; i<nset; i++)
  {
    tbeg = DBL_MAX;
    tend = 0.0;
    for(j=0; j<dataset[i].nlset; j++)
    {
      tbeg = fmin(tbeg, dataset[i].line[j].t[0]);
      tend = fmax(tend, dataset[i].line[j].t[dataset[i].line[j].n-1]);
    }
    tspan = tend - tbeg;
    tbeg -= 0.05*tspan;
    tend += 0.05*tspan;
    tspan = tend - tbeg;

    /* time nodes of continuum */
    for(j=0; j<nall[i][0]; j++)
      tall[i][j] = tspan/(nall[i][0]-1.0) * j + tbeg;

    /* time nodes of lines */
    np = nall[i][0];
    for(j=0; j<dataset[i].nlset; j++)
    {
      for(k=0; k<nall[i][1+j]; k++)
      {
        tall[i][np+k] = tspan/(nall[i][1+j]-1.0) * k + tbeg;
      }
      np += nall[i][1+j];
    }
  }

  for(i=0; i<nset; i++)
  { 
    for(k=0; k<nall[i][0]; k++)
    {
      fall_best[i][k] = 0.0;
      fall_std[i][k] = 0.0;
      feall_max[i][k] = 0.0;
    }
      
    np = nall[i][0];
    for(j=0; j<dataset[i].nlset; j++)
    {
      for(k=0; k<nall[i][1+j]; k++)
      {
        fall_best[i][np+k] = 0.0;
        fall_std[i][np+k] = 0.0;
        feall_max[i][np+k] = 0.0;
      }          
      np += nall[i][1+j];
    } 

    for(k=0; k<(1+nlset_max)*nq; k++)
    {
      yq_best[i][k] = 0.0;
      yq_std[i][k] = 0.0;
    } 
  }
  
  post_model = (double *)posterior_sample_task;
  for(m=0; m<num_ps_task; m++)
  {
    printf("# sample %d on task %d\n", m, thistask);
    
    for(i=0; i<nset; i++)
    {
      /* reconstuct all the light curves */
      reconstruct_line_from_varmodel_vmap((void *)post_model, i, nall[i], tall[i], fall[i], feall[i], yq); 
      
      for(k=0; k<nall[i][0]; k++)
      {
        fall_best[i][k] += fall[i][k];
        fall_std[i][k] += fall[i][k]*fall[i][k];
        feall_max[i][k] = fmax(feall_max[i][k], feall[i][k]);
      }
      /* reconstructed lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
        {
          fall_best[i][np+k] += fall[i][np+k];
          fall_std[i][np+k] += fall[i][np+k]*fall[i][np+k];
          feall_max[i][np+k] = fmax(feall_max[i][np+k], feall[i][np+k]);
        }          
        np += nall[i][1+j];
      }  
      
      for(k=0; k<nq*(1+dataset[i].nlset); k++)
      {
        yq_best[i][k] += yq[k];
        yq_std[i][k] += yq[k]*yq[k];
      }
    }

    post_model += num_params;
  }

  /* now gather the reconstruction */
  for(i=0; i<nset; i++)
  {
    MPI_Reduce(fall_best[i], fall_best_buf[i], ntall[i], MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
    MPI_Reduce(fall_std[i], fall_std_buf[i], ntall[i], MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
    MPI_Reduce(feall_max[i], feall_max_buf[i], ntall[i], MPI_DOUBLE, MPI_MAX, roottask, MPI_COMM_WORLD);
    MPI_Reduce(yq_best[i], yq_best_buf[i], nq*(1+dataset[i].nlset), MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
    MPI_Reduce(yq_std[i], yq_std_buf[i], nq*(1+dataset[i].nlset), MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
  }

  if(thistask == roottask)
  {
    /* write headers */
    for(i=0; i<nset; i++)
    {
      fprintf(fp, "# %d",nall[i][0]);
      for(j=0;  j < dataset[i].nlset; j++)
        fprintf(fp, ":%d", nall[i][1+j]);
      fprintf(fp, "\n");
    }

    for(i=0; i<nset; i++)
    {
      for(k=0; k<ntall[i]; k++)
      {
        fall_best_buf[i][k] /= num_ps;
        fall_std_buf[i][k] /= num_ps;
        fall_std_buf[i][k] = sqrt(fall_std_buf[i][k] - fall_best_buf[i][k]*fall_best_buf[i][k]);
        fall_std_buf[i][k] = fmax(fall_std_buf[i][k], feall_max_buf[i][k]);
      }

      /* output reconstructed continuum */
      for(k=0; k<nall[i][0]; k++)
        fprintf(fp, "%e %e %e\n", tall[i][k], fall_best_buf[i][k] * dataset[i].con.scale, fall_std_buf[i][k] * dataset[i].con.scale);
      fprintf(fp, "\n");

      /* output reconstructed lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
          fprintf(fp, "%e %e %e\n", tall[i][np+k], fall_best_buf[i][np+k] * dataset[i].line[j].scale, fall_std_buf[i][np+k] * dataset[i].line[j].scale);
        fprintf(fp, "\n");
        np += nall[i][1+j];
      }  

      /* output long-term trend yq */
      for(j=0; j<(1+dataset[i].nlset)*nq; j++)
      {
        yq_best_buf[i][j] /= num_ps;
        yq_std_buf[i][j] /= num_ps;
        yq_std_buf[i][j] = sqrt(yq_std_buf[i][j] - yq_best_buf[i][j]*yq_best_buf[i][j]);
      }
      for(j=0; j<1*nq; j++)
        fprintf(fpq, "%e %e\n", yq_best_buf[i][j]* dataset[i].con.scale, yq_std_buf[i][j]* dataset[i].con.scale);
      for(j=0; j<dataset[i].nlset; j++)
        for(k=0; k<nq; k++)
          fprintf(fpq, "%e %e\n", yq_best_buf[i][(1+j)*nq+k]* dataset[i].line[j].scale, yq_std_buf[i][(1+j)*nq+k]* dataset[i].line[j].scale);
      fprintf(fpq, "\n");
      
      if(parset.flag_trend > 0)
      {
        printf("Longterm q of dataset %d: Val and Err\n", i);
        for(j=0; j<1*nq; j++)
        {
          printf("%e %e\n", yq_best_buf[i][j]*dataset[i].con.scale, yq_std_buf[i][j]* dataset[i].con.scale);
        }
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nq; k++)
            printf("%e %e\n", yq_best_buf[i][(1+j)*nq+k] * dataset[i].line[j].scale, yq_std_buf[i][(1+j)*nq+k]* dataset[i].line[j].scale);
        }
      }
    }

    fclose(fp);
    fclose(fpq);
  }

  free(posterior_sample_task);
  if(thistask == roottask)
  {
    free(posterior_sample_all);
  }

  for(i=0; i<nset; i++)
  {
    free(tall[i]);
    free(fall[i]);
    free(feall[i]);
    free(nall[i]);
    free(feall_max[i]);
    free(fall_best[i]);
    free(fall_std[i]);

    free(yq_best[i]);
    free(yq_std[i]);

    free(feall_max_buf[i]);
    free(fall_best_buf[i]);
    free(fall_std_buf[i]);

    free(yq_best_buf[i]);
    free(yq_std_buf[i]);
  }
  free(tall);
  free(fall);
  free(feall);
  free(feall_max);
  free(fall_best);
  free(fall_std);
  free(nall);
  free(ntall);
  free(yq);
  free(yq_best);
  free(yq_std);

  free(feall_max_buf);
  free(fall_best_buf);
  free(fall_std_buf);

  free(yq_best_buf);
  free(yq_std_buf);

  return;
}

/* 
 * reconstruct line lcs acording to the time span of the continuum.
 */
void output_reconstruction_vmap()
{
  if(thistask == roottask)
  {
    printf("start reconstruction...\n");

    FILE *fp, *fp_sample, *fpq;
    char fname[200];
    int i, j, k, m;
    double **tall, **fall, **feall, **feall_max, **fall_best, **fall_std, *yq, **yq_best, **yq_std;
    int **nall, *ntall, np;
    double tspan;
    
    int num_ps, size_of_modeltype;
    void *post_model;
    char posterior_sample_file[MICA_MAX_STR_LENGTH];

    /* time span of reconstruction */
    double tbeg, tend;

    size_of_modeltype = num_params * sizeof(double);

    /* get file name of posterior sample file */
    dnest_get_posterior_sample_file(posterior_sample_file);
    /* open file for posterior sample */
    fp_sample = fopen(posterior_sample_file, "r");
    if(fp_sample == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    /* read number of points in posterior sample */
    if(fscanf(fp_sample, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }

    sprintf(fname, "%s/%s%s", parset.file_dir, "data/pall.txt", postfix);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    fprintf(fp, "# %d\n", nset);

    sprintf(fname, "%s/%s%s", parset.file_dir, "data/trend.txt", postfix);
    fpq = fopen(fname, "w");
    if(fpq == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
 
    nall = malloc(nset*sizeof(int *));
    tall = malloc(nset*sizeof(double *));
    fall = malloc(nset*sizeof(double *));
    feall = malloc(nset*sizeof(double *));
    feall_max = malloc(nset*sizeof(double *));
    fall_best = malloc(nset*sizeof(double *));
    fall_std = malloc(nset*sizeof(double *));
    ntall = malloc(nset*sizeof(int));
    yq_best = malloc(nset*sizeof(double *));
    yq_std = malloc(nset*sizeof(double *));
    
    /* first set the number of reconstruction points */
    for(i=0; i<nset; i++)
    {
      nall[i] = malloc((1+nlset_max) * sizeof(int));
      
      nall[i][0] = 0;
      for(j=0; j<dataset[i].nlset; j++)
      {
        /* number < 500 */
        nall[i][0] = fmax(nall[i][0], fmin(dataset[i].line[j].n*nscale, 500));
      }

      for(j=0; j<dataset[i].nlset; j++)
        nall[i][1+j] = nall[i][0];
      
      /* compute total number of points */
      ntall[i] = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
        ntall[i] += nall[i][1+j];
    }

    for(i=0; i<nset; i++)
    { 
      tall[i] = malloc(ntall[i] * sizeof(double));
      fall[i] = malloc(ntall[i] * sizeof(double));
      feall[i] = malloc(ntall[i] * sizeof(double));
      feall_max[i] = malloc(ntall[i] * sizeof(double));
      fall_best[i] = malloc(ntall[i] * sizeof(double));
      fall_std[i] = malloc(ntall[i] * sizeof(double));

      yq_best[i] = malloc(nq*(nlset_max)*sizeof(double));
      yq_std[i] = malloc(nq*(nlset_max)*sizeof(double));
    }
    
    post_model = malloc(size_of_modeltype);
    yq = malloc(nq*(nlset_max)*sizeof(double));

    for(i=0; i<nset; i++)
    {
      tbeg = DBL_MAX;
      tend = 0.0;
      for(j=0; j<dataset[i].nlset; j++)
      {
        tbeg = fmin(tbeg, dataset[i].line[j].t[0]);
        tend = fmax(tend, dataset[i].line[j].t[dataset[i].line[j].n-1]);
      }
      tspan = tend - tbeg;
      tbeg -= 0.05*tspan;
      tend += 0.05*tspan;
      tspan = tend - tbeg;

      /* time nodes of continuum */
      for(j=0; j<nall[i][0]; j++)
        tall[i][j] = tspan/(nall[i][0]-1.0) * j + tbeg;

      /* time nodes of lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
        {
          tall[i][np+k] = tspan/(nall[i][1+j]-1.0) * k + tbeg;
        }
        np += nall[i][1+j];
      }
    }

    for(i=0; i<nset; i++)
    { 
      for(k=0; k<nall[i][0]; k++)
      {
        fall_best[i][k] = 0.0;
        fall_std[i][k] = 0.0;
        feall_max[i][k] = 0.0;
      }
        
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
        {
          fall_best[i][np+k] = 0.0;
          fall_std[i][np+k] = 0.0;
          feall_max[i][np+k] = 0.0;
        }          
        np += nall[i][1+j];
      } 

      for(k=0; k<(nlset_max)*nq; k++)
      {
        yq_best[i][k] = 0.0;
        yq_std[i][k] = 0.0;
      } 
    }

    for(m=0; m<num_ps; m++)
    {
      printf("# sample %d\n", m);
      
      // read sample
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp_sample, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp_sample, "\n");
      
      for(i=0; i<nset; i++)
      {
        /* reconstuct all the light curves */
        reconstruct_line_from_varmodel_vmap(post_model, i, nall[i], tall[i], fall[i], feall[i], yq); 

        for(k=0; k<nall[i][0]; k++)
        {
          fall_best[i][k] += fall[i][k];
          fall_std[i][k] += fall[i][k]*fall[i][k];
          feall_max[i][k] = fmax(feall_max[i][k], feall[i][k]);
        }
        /* reconstructed lines */
        np = nall[i][0];
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nall[i][1+j]; k++)
          {
            fall_best[i][np+k] += fall[i][np+k];
            fall_std[i][np+k] += fall[i][np+k]*fall[i][np+k];
            feall_max[i][np+k] = fmax(feall_max[i][np+k], feall[i][np+k]);
          }          
          np += nall[i][1+j];
        }  
        
        for(k=0; k<nq*(dataset[i].nlset); k++)
        {
          yq_best[i][k] += yq[k];
          yq_std[i][k] += yq[k]*yq[k];
        } 
      }
    }
    
    /* write headers */
    for(i=0; i<nset; i++)
    {
      fprintf(fp, "# %d",nall[i][0]);
      for(j=0;  j < dataset[i].nlset; j++)
        fprintf(fp, ":%d", nall[i][1+j]);
      fprintf(fp, "\n");
    }

    for(i=0; i<nset; i++)
    { 
      for(k=0; k<nall[i][0]; k++)
      {
        fall_best[i][k] /= num_ps;
        fall_std[i][k] /= num_ps;
        fall_std[i][k] = sqrt(fall_std[i][k] - fall_best[i][k]*fall_best[i][k]);
        fall_std[i][k] = fmax(fall_std[i][k], feall_max[i][k]);
      }
        
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
        {
          fall_best[i][np+k] /= num_ps;
          fall_std[i][np+k] /= num_ps;
          fall_std[i][np+k] = sqrt(fall_std[i][np+k] - fall_best[i][np+k]*fall_best[i][np+k]);
          fall_std[i][np+k] = fmax(fall_std[i][np+k], feall_max[i][np+k]);
        }          
        np += nall[i][1+j];
      }
      
      /* output reconstructed continuum */
      for(k=0; k<nall[i][0]; k++)
        fprintf(fp, "%e %e %e\n", tall[i][k], fall_best[i][k] * dataset[i].con.scale, fall_std[i][k] * dataset[i].con.scale);
      fprintf(fp, "\n");

      /* output reconstructed lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
          fprintf(fp, "%e %e %e\n", tall[i][np+k], fall_best[i][np+k] * dataset[i].line[j].scale, fall_std[i][np+k] * dataset[i].line[j].scale);
        fprintf(fp, "\n");
        np += nall[i][1+j];
      }  

      /* output long-term trend yq */
      for(j=0; j<(dataset[i].nlset)*nq; j++)
      {
        yq_best[i][j] /= num_ps;
        yq_std[i][j] /= num_ps;
        yq_std[i][j] = sqrt(yq_std[i][j] - yq_best[i][j]*yq_best[i][j]);
      }
      for(j=0; j<1*nq; j++)
        fprintf(fpq, "%e %e\n", 0.0*dataset[i].con.scale, 0.0*dataset[i].con.scale);
      for(j=0; j<dataset[i].nlset; j++)
        for(k=0; k<nq; k++)
          fprintf(fpq, "%e %e\n", yq_best[i][(j)*nq+k]* dataset[i].line[j].scale, yq_std[i][(j)*nq+k]* dataset[i].line[j].scale);
      fprintf(fpq, "\n");
      
      if(parset.flag_trend > 0)
      {
        printf("Longterm q of dataset %d: Val and Err\n", i);
        for(j=0; j<1*nq; j++)
        {
          printf("%e %e\n", 0.0*dataset[i].con.scale, 0.0*dataset[i].con.scale);
        }
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nq; k++)
            printf("%e %e\n", yq_best[i][(j)*nq+k] * dataset[i].line[j].scale, yq_std[i][(j)*nq+k]* dataset[i].line[j].scale);
        }
      }
    }

    fclose(fp);
    fclose(fp_sample);
    fclose(fpq);

    for(i=0; i<nset; i++)
    {
      free(tall[i]);
      free(fall[i]);
      free(feall[i]);
      free(nall[i]);
      free(feall_max[i]);
      free(fall_best[i]);
      free(fall_std[i]);

      free(yq_best[i]);
      free(yq_std[i]);
    }
    free(tall);
    free(fall);
    free(feall);
    free(feall_max);
    free(fall_best);
    free(fall_std);
    free(nall);
    free(ntall);
    free(post_model);
    free(yq);
    free(yq_best);
    free(yq_std);

    printf("end reconstruction...\n");
  }

  return;
}

/*!
 * matrix operation A^-1 x B is implemented by calling function 
 *    multiply_mat_MN_inverseA()
 * 
 */
void reconstruct_line_from_varmodel_vmap(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq, *yave;
  int i, j, k, m, info, idx, *ipiv;
  double *PEmat1, *PEmat2, *PEmat3, *PEmat4;
  int nall_data, nqall, ntall, np, ntall_max;
  double *fall_data;
  double sigma, tau, *pm=(double *)model;

  ipiv = workspace_ipiv;

  idx = idx_con_pm[nds];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]); //* sqrt(tau);

  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  nqall = nq * (dataset[nds].nlset);
  nall_data = alldata[nds].n;
  fall_data = alldata[nds].f;

  Larr = workspace;
  ybuf = Larr + nall_data * nqall;
  y = ybuf + nall_data * nqall;
  Cq = y + nall_data;
  yq = Cq + nqall*nqall;
  yave = yq + nqall;
  yuq = yave + nall_data; 
  Larr_rec = yuq + ntall;
  
  /* determine the maximum size */
  ntall_max = ntall>nall_max?ntall:nall_max;
  
  PEmat1 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat2 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat3 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat4 = malloc(ntall_max * ntall_max * sizeof(double));

  set_covar_Pmat_data_line_vmap(model, nds);
  
  set_covar_Umat_line_vmap(model, nds, nall, tall);
  
  set_covar_Amat_line_vmap(model, nds, nall, tall);
  
  /* no continuum data point */
  np = 0;
  for(j=0; j<dataset[nds].nlset; j++)
  {
    for(m=0; m<dataset[nds].line[j].n; m++)
    {
      for(i=0; i<nqall; i++)
        Larr[(np+m)*nqall + i ]  = 0.0;
        
      Larr[(np+m)*nqall + j*nq + 0] = 1.0;
      for(i=1; i<nq; i++)
        Larr[(np+m)*nqall + j*nq + i] = pow(dataset[nds].line[j].t[m], i);
    }
    np += dataset[nds].line[j].n;
  }

  memcpy(PEmat1, PCmat, nall_data*nall_data*sizeof(double));
  memcpy(PEmat2, Larr, nall_data*nqall*sizeof(double));

  multiply_mat_MN_inverseA(PEmat1, PEmat2, nall_data, nqall, ipiv); // Tmat2 = C^-1 * L;  NxNq

  multiply_mat_MN_transposeA(Larr, PEmat2, Cq, nqall, nqall, nall_data); // ICq = L^T*C^-1*L; NqxNq
  multiply_mat_MN_transposeA(PEmat2, fall_data, yq, nqall, 1, nall_data); // yq = L^T*C^-1*y;  Nqx1
  memcpy(PEmat1, Cq, nqall*nqall*sizeof(double));
  multiply_mat_MN_inverseA(PEmat1, yq, nqall, 1, ipiv); // yq = (L^T*C^-1*L)^-1 * L^T*C^-1*y; Nqx1
  memcpy(yqall, yq, nqall*sizeof(double));

  multiply_mat_MN(Larr, yq, yave, nall_data, 1, nqall); // yave = L * q; Nx1
  for(i=0; i<nall_data; i++)
    y[i] = fall_data[i] - yave[i];

  memcpy(PEmat1, PCmat, nall_data*nall_data*sizeof(double));
  memcpy(PEmat2, y, nall_data*sizeof(double));
  multiply_mat_MN_inverseA(PEmat1, PEmat2, nall_data, 1, ipiv);  // C^-1 * (y - Lq)

  multiply_matvec_MN(USmat, ntall, nall_data, PEmat2, fall); // S * C^-1 * (y - Lq)

  /* note that no continuum data point, Larr_rec is set to be zero */
  for(i=0;i<nall[0];i++)
  {
    for(j=0; j<nqall; j++)
      Larr_rec[i*nqall + j] = 0.0;
    
    Larr_rec[i*nqall + 0] = 0.0;
  }
  
  /* now for lines */
  np = nall[0];
  for(k=0; k<dataset[nds].nlset; k++)
  {
    for(i=0; i<nall[1+k]; i++)
    {
      for(j=0; j<nqall; j++)
      {
        Larr_rec[(np+i)*nqall + j] = 0.0;
      }

      Larr_rec[(np+i)*nqall + k*nq + 0] = 1.0;
      for(j=1; j<nq; j++)
        Larr_rec[(np+i)*nqall + k*nq + j] = pow(tall[np+i], j);
    }
    np += nall[1+k];
  }

  multiply_matvec_MN(Larr_rec, ntall, nqall, yq, yuq);

  for(i=0; i<ntall; i++)
    fall[i] += yuq[i];

  /* Transpose of USmat */
  for(i=0; i<ntall; i++)
    for(j=0; j<nall_data; j++)
      USmatT[j*ntall + i] = USmat[i*nall_data + j];

  memcpy(PEmat1, PCmat, nall_data*nall_data*sizeof(double));
  memcpy(PEmat2, USmatT, ntall*nall_data*sizeof(double));

  multiply_mat_MN_inverseA(PEmat1, PEmat2, nall_data, ntall, ipiv); // C^-1 x S; NdxN
  multiply_mat_MN(USmat, PEmat2, PEmat1, ntall, ntall, nall_data); // S x C^-1 x S; NxN

  /* S x C^-1 x L; NxNq */
  multiply_mat_MN_transposeA(PEmat2, Larr, PEmat3, ntall, nqall, nall_data); 
  /* S x C^-1 x L - L; NxNq */
  for(i=0; i<ntall*nqall; i++)PEmat3[i] -= Larr_rec[i];
  
  inverse_mat(Cq, nqall, &info, ipiv);

  multiply_mat_MN(PEmat3, Cq, PEmat2, ntall, nqall, nqall);
  /* (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  multiply_mat_MN_transposeB(PEmat2, PEmat3, PEmat4, ntall, ntall, nqall);

  /* continuum, no q involved */  
  for(i=0; i<nall[0]; i++)
  {
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat1[i*ntall+i]);
  }
  for(i=nall[0]; i<ntall; i++)
  {
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat1[i*ntall+i] + PEmat4[i*ntall+i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  return;
}

/*!
 * this function sets the covariance matrix at data time points 
 *
 * k is the index of dataset
 */
void set_covar_Pmat_data_line_vmap(const void *model, int k)
{
  double t1, t2;
  double sigma, tau, syserr, syserr_line, error;
  int i, j, m, l, nall, np, npline, idx;
  double *pm = (double *)model;
  
  /* total number of point of this dataset */
  nall = alldata[k].n;

  /* set variability parameters */
  idx = idx_con_pm[k];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]); // * sqrt(tau); 

  /* no continuum data points */

  /* line - line */
  np = 0;
  for(j=0; j<dataset[k].nlset; j++)
  {

    idx = idx_line_pm[k][j];
    syserr_line = (exp(pm[idx]) - 1.0) * dataset[k].line[j].error_mean;

    for(i=0; i<dataset[k].line[j].n; i++)
    {
      /* line self */
      t1 = dataset[k].line[j].t[i];
      for(m=0; m<i; m++)
      {
        t2 = dataset[k].line[j].t[m];
        
        PCmat[(np + i)*nall+ (np + m) ] = PCmat[(np+m)*nall+(np+i)] =
        PSmat[(np + i)*nall+ (np + m) ] = PSmat[(np+m)*nall+(np+i)] = Sll(t1, t2, model, k, j);
      }
      PSmat[(np + i)*nall+(np +i)] = Sll(t1, t1, model, k, j);
      error = dataset[k].line[j].fe[i] * dataset[k].line[j].fe[i] + syserr_line*syserr_line;
      PCmat[(np + i)*nall+(np +i)] = PSmat[(np + i)*nall+(np +i)] + error/sigma/sigma;

      /* between lines */
      npline = np + dataset[k].line[j].n;
      for(l=j+1; l<dataset[k].nlset; l++)
      {
        for(m = 0; m < dataset[k].line[l].n; m++)
        {
          t2 = dataset[k].line[l].t[m];

          PCmat[(np+i)*nall + npline + m ] = PCmat[ (npline + m) * nall + np + i ] =
          PSmat[(np+i)*nall + npline + m ] = PSmat[ (npline + m) * nall + np + i ] = Sll2(t1, t2, model, k, j, l);
        }

        npline += dataset[k].line[l].n;
      }
    }
    np += dataset[k].line[j].n;
  }

  return;
}

/*!
 * this function sets the covariance matrix at time of data points and reconstruction points
 */
void set_covar_Umat_line_vmap(const void *model, int nds, int *nall, double *tall)
{
  double taud, t1, t2;
  int i, j, k, m, ntall, nall_data, np, npline, idx;
  double *pm = (double *)model;

  nall_data = alldata[nds].n;
  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  /* continuum - continuum/line */
  for(i=0; i<nall[0]; i++)
  {
    /* reconstruction point */
    t1 = tall[i];
   
    /* no continuum data points */

    /* lines */
    np = 0;
    for(k=0; k<dataset[nds].nlset; k++)
    {
      for(j=0; j < dataset[nds].line[k].n; j++)
      {
        t2 = dataset[nds].line[k].t[j];
        USmat[ i * nall_data + j + np ] = Slc(t1, t2, model, nds, k);
      }
      np += dataset[nds].line[k].n;
    }
  }
  
  /* line - line */
  npline = nall[0];
  for(k = 0; k<dataset[nds].nlset; k++)
  {
    for(i=0; i<nall[1+k]; i++)
    {
      t1 = tall[npline + i];

      /* line different */
      np = 0;
      for(m=0; m<k; m++)
      {
        for(j=0; j < dataset[nds].line[m].n; j++)
        {
          t2 = dataset[nds].line[m].t[j];
          USmat[ (i + npline) * nall_data + j + np ] = Sll2(t1, t2, model, nds, k, m);
        }
        np += dataset[nds].line[m].n;
      }
      
      /* line self */
      for(j=0; j<dataset[nds].line[k].n; j++)
      {
        t2 = dataset[nds].line[k].t[j];
        USmat[ (i+npline) * nall_data + j + np ] = Sll(t1, t2, model, nds, k);
      }
      np += dataset[nds].line[k].n;

      /* line different */
      for(m=k+1; m<dataset[nds].nlset; m++)
      {
        for(j=0; j < dataset[nds].line[m].n; j++)
        {
          t2 = dataset[nds].line[m].t[j];
          USmat[ (i + npline) * nall_data + j + np ] = Sll2(t1, t2, model, nds, k, m);
        }
        np += dataset[nds].line[m].n;
      }
    }
    npline += nall[1+k];
  }
  return;
}

/*!
 * this function sets the covariance matrix at time of reconstruction points
 *
 * exactly same with set_covar_Amat_line()
 */
void set_covar_Amat_line_vmap(const void *model, int nds, int *nall, double *tall)
{
  double t1, t2, taud;
  int i, j, k, m, ntall, np, npline, idx;
  double *pm=(double *)model;
  
  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  /* continuum - continuum/lines */
  for(i=0; i<nall[0]; i++)
  {
    t1 = tall[i];

    /* continuum */
    for(j=0; j<i; j++)
    {
      t2 = tall[j];
      ASmat[i * ntall + j] = ASmat[j * ntall + i] =  exp( - fabs(t1-t2)/taud);
    }
    ASmat[i * ntall + i] = 1.0;
    
    /* lines */
    np = nall[0];
    for(k=0; k<dataset[nds].nlset; k++)
    {
      for(j=0; j<nall[1+k]; j++)
      {
        t2 = tall[np + j];
        ASmat[i * ntall + np+j] = ASmat[(np+j)*ntall + i] = Slc(t1, t2, model, nds, k);
      }
      np += nall[1+k];
    }
  }

  /* line - line */
  np = nall[0];
  for(k=0; k<dataset[nds].nlset; k++)
  {
    for(j = 0; j<nall[1+k]; j++)
    {
      t1 = tall[np + j];

      /* line self */
      for(i=0; i<j; i++)
      {
        t2 = tall[np + i];
        ASmat[(j+np)*ntall + np+i] = ASmat[(i+np)*ntall + np+j] = Sll(t1, t2, model, nds, k);
      }
      ASmat[(j+np)*ntall + (j+np)] = Sll(t1, t1, model, nds, k);
      
      /* line different */
      npline = np + nall[1+k];
      for(m = k+1; m<dataset[nds].nlset; m++)
      {
        for(i=0; i<nall[1+m]; i++)
        {
          t2 = tall[npline + i];
          ASmat[(j+np)*ntall + npline + i] = ASmat[ (npline + i)*ntall + j+np ] = Sll2(t1, t2, model, nds, k, m) ;
        }
        npline += nall[1+m];
      }
    }
    np += nall[1+k];
  }
  return;
}

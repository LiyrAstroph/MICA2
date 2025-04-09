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

#include "dnest_line.h"
#include "allvars.h"
#include "dnest.h"

#include "proto.h"

void mc_line()
{
  int i, j, argc=0, jzmax=0, narg, narg0;
  double logz_max;
  char **argv;

  argv = malloc(11*sizeof(char *));
  for(i=0; i<11; i++)
  {
    argv[i] = malloc(MICA_MAX_STR_LENGTH*sizeof(char));
  }
  /* setup argc and argv 
   * 0: dnest
   * 1: -s
   * 2: file name 
   */
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  //strcpy(argv[argc], parset.file_dir);
  //strcat(argv[argc++], "/data/restart_dnest1d.txt");
  argc++;
  narg=3;
  

  /* level-dependent sampling */
  {
    strcpy(argv[argc++], "-l");
    narg++;
  }

  if(flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
    narg++;
  }
  
  // sample tag 
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "1d");
  narg+=2;

  mc_line_init();
  
  narg0 = narg;
  logz_max = -DBL_MAX;
  for(j=0; j<parset.num_gaussian_diff; j++)
  {
    num_gaussian = parset.num_gaussian_low + j;
    
    type_lag_prior_pr = 1; 
    if(parset.type_lag_prior == 0 && num_gaussian > 1)
    {
      type_lag_prior_pr = 0;
    }
    
    if(thistask == roottask)
    {
      printf("# number of Gaussian/Tophat: %d\n", num_gaussian);
    }
    
    sprintf(postfix, "_%d", num_gaussian);
    strcpy(argv[argc], "-x");
    strcpy(argv[argc+1], postfix);
    narg=narg0+2;
    
    /* restart save file name */
    strcpy(argv[2], parset.file_dir);
    strcat(argv[2], "/data/restart_dnest1d.txt");
    strcat(argv[2], postfix);

    if(flag_restart == 1)
    {
      /* restart file name */
      strcpy(argv[argc+2], "-r");
      strcpy(argv[argc+3], parset.file_dir);
      strcat(argv[argc+3], "/");
      strcat(argv[argc+3], "data/restart_dnest1d.txt");
      strcat(argv[argc+3], postfix);
      narg+=2;
    }

    logz_arr[j] = dnest_line(narg, argv);
    if(logz_max < logz_arr[j])
    {
      logz_max = logz_arr[j];
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

  mc_line_end();

  /* clear up argv */
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

/*
 * do reconstrunction in parallel
 */
void output_reconstruction_parallel()
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
    nall[i] = (int *)malloc((1+nlset_max) * sizeof(int));

    nall[i][0] = (int) fmin(dataset[i].con.n*nscale, parset.nd_rec);
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
    tbeg = dataset[i].con.t[0];
    tend = dataset[i].con.t[dataset[i].con.n-1];
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
    if(parset.model == pmap)
    {
      transform_response_ratio_inplace((void *)post_model);
    }
    
    for(i=0; i<nset; i++)
    {
      /* reconstuct all the light curves */
      reconstruct_line_from_varmodel3((void *)post_model, i, nall[i], tall[i], fall[i], feall[i], yq); 
      
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
void output_reconstruction()
{
  if(thistask == roottask)
  {
    printf("start reconstruction...\n");

    FILE *fp, *fp_sample, *fpq, *fp_each;
    char fname[200];
    int i, j, k, m;
    double **tall, **fall, **feall, **feall_max, **fall_best, **fall_std, *yq, **yq_best, **yq_std;
    int **nall, *ntall, np;
    double tspan;
    
    int num_ps, size_of_modeltype;
    void *post_model, *post_model_trans; 
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

    sprintf(fname, "%s/%s%s", parset.file_dir, "data/pall_each.txt", postfix);
    fp_each = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    fprintf(fp_each, "# %d\n", nset);

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

      nall[i][0] = (int) fmin(dataset[i].con.n*nscale, parset.nd_rec);
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

      yq_best[i] = malloc(nq*(1+nlset_max)*sizeof(double));
      yq_std[i] = malloc(nq*(1+nlset_max)*sizeof(double));
    }
    
    post_model = malloc(size_of_modeltype);
    yq = malloc(nq*(1+nlset_max)*sizeof(double));

    /* in pmap, response ratio is used, need to transform into normal values as in gmodel */
    if(parset.model == pmap)
    {
      post_model_trans = malloc(size_of_modeltype);
    }
    else
    {
      post_model_trans = post_model;
    }

    for(i=0; i<nset; i++)
    {
      tbeg = dataset[i].con.t[0];
      tend = dataset[i].con.t[dataset[i].con.n-1];
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
      
      if(parset.model == pmap)
      {
        transform_response_ratio(post_model, post_model_trans);
      }

      for(i=0; i<nset; i++)
      {
        /* reconstuct all the light curves */
        reconstruct_line_from_varmodel3(post_model_trans, i, nall[i], tall[i], fall[i], feall[i], yq); 
        
        /* output reconstructed lines */
        np = nall[i][0];
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nall[i][1+j]; k++)
            fprintf(fp_each, "%e %e %e\n", tall[i][np+k], fall[i][np+k] * dataset[i].line[j].scale, feall[i][np+k] * dataset[i].line[j].scale);
          fprintf(fp_each, "\n");
          np += nall[i][1+j];
        }  

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
      for(j=0; j<(1+dataset[i].nlset)*nq; j++)
      {
        yq_best[i][j] /= num_ps;
        yq_std[i][j] /= num_ps;
        yq_std[i][j] = sqrt(yq_std[i][j] - yq_best[i][j]*yq_best[i][j]);
      }
      for(j=0; j<1*nq; j++)
        fprintf(fpq, "%e %e\n", yq_best[i][j]* dataset[i].con.scale, yq_std[i][j]* dataset[i].con.scale);
      for(j=0; j<dataset[i].nlset; j++)
        for(k=0; k<nq; k++)
          fprintf(fpq, "%e %e\n", yq_best[i][(1+j)*nq+k]* dataset[i].line[j].scale, yq_std[i][(1+j)*nq+k]* dataset[i].line[j].scale);
      fprintf(fpq, "\n");
      
      if(parset.flag_trend > 0)
      {
        printf("Longterm q of dataset %d: Val and Err\n", i);
        for(j=0; j<1*nq; j++)
        {
          printf("%e %e\n", yq_best[i][j]*dataset[i].con.scale, yq_std[i][j]* dataset[i].con.scale);
        }
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nq; k++)
            printf("%e %e\n", yq_best[i][(1+j)*nq+k] * dataset[i].line[j].scale, yq_std[i][(1+j)*nq+k]* dataset[i].line[j].scale);
        }
      }
    }

    fclose(fp);
    fclose(fp_sample);
    fclose(fpq);
    fclose(fp_each);

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
    if(parset.model == pmap)
    {
      free(post_model_trans);
    }
    free(yq);
    free(yq_best);
    free(yq_std);

    printf("end reconstruction...\n");
  }

  return;
}

/* 
 * reconstruct line lcs acording to their own time spans.
 */
void output_reconstruction2()
{
  if(thistask == roottask)
  {
    FILE *fp, *fp_sample, *fpq;
    char fname[200];
    int i, j, k, m;
    double **tall, **fall, **feall, **feall_max, **fall_best, **fall_std, *yq, **yq_best, **yq_std;
    int **nall, *ntall, np;
    double tspan;
    
    int num_ps, size_of_modeltype;
    void *post_model, *post_model_trans;
    char posterior_sample_file[MICA_MAX_STR_LENGTH];

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
    for(i=0; i<nset; i++)
    {
      nall[i] = malloc((1+nlset_max) * sizeof(int));
      tall[i] = malloc(nall_max*nscale * sizeof(double));
      fall[i] = malloc(nall_max*nscale * sizeof(double));
      feall[i] = malloc(nall_max*nscale * sizeof(double));
      feall_max[i] = malloc(nall_max*nscale * sizeof(double));
      fall_best[i] = malloc(nall_max*nscale * sizeof(double));
      fall_std[i] = malloc(nall_max*nscale * sizeof(double));

      yq_best[i] = malloc(nq*(1+nlset_max)*sizeof(double));
      yq_std[i] = malloc(nq*(1+nlset_max)*sizeof(double));
    }
    
    post_model = malloc(size_of_modeltype);
    yq = malloc(nq*(1+nlset_max)*sizeof(double));
    
    /* in pmap, response ratio is used, need to transform into normal values as in gmodel */
    if(parset.model == pmap)
    {
      post_model_trans = malloc(size_of_modeltype);
    }
    else
    {
      post_model_trans = post_model;
    }

    for(i=0; i<nset; i++)
    {
      /* time nodes of continuum */
      nall[i][0] = (int)fmin(dataset[i].con.n * nscale, parset.nd_rec);
      tspan = dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0];
      for(j=0; j<nall[i][0]; j++)
        tall[i][j] = (tspan+0.1*tspan)/(nall[i][0]-1.0) * j 
                        + dataset[i].con.t[0]-0.05*tspan;

      /* time nodes of lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        nall[i][1+j] = (int)fmin(dataset[i].line[j].n * nscale, parset.nd_rec);
        tspan = dataset[i].line[j].t[dataset[i].line[j].n-1] - dataset[i].line[j].t[0];
        for(k=0; k<nall[i][1+j]; k++)
        {
          tall[i][np+k] = (tspan+0.1*tspan)/(nall[i][1+j]-1.0) * k 
                         + dataset[i].line[j].t[0]-0.05*tspan;
        }
        np += nall[i][1+j];
      }

      /* compute total number of points */
      ntall[i] = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
        ntall[i] += nall[i][j];
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

    for(m=0; m<num_ps; m++)
    {
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
      
      if(parset.model == pmap)
      {
        transform_response_ratio(post_model, post_model_trans);
      }

      for(i=0; i<nset; i++)
      {
        /* reconstuct all the light curves */
        reconstruct_line_from_varmodel3(post_model_trans, i, nall[i], tall[i], fall[i], feall[i], yq); 

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
    }

    /* first output number of points of reconstructions */
    for(i=0; i<nset; i++)
    {
      fprintf(fp, "# %d",nall[i][0]);
      for(j=0;  j < dataset[i].nlset; j++)
        fprintf(fp, ":%d", nall[i][1+j]);
      fprintf(fp, "\n");
    }
    /* then output reconstructions */
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
      for(j=0; j<(1+dataset[i].nlset)*nq; j++)
      {
        yq_best[i][j] /= num_ps;
        yq_std[i][j] /= num_ps;
        yq_std[i][j] = sqrt(yq_std[i][j] - yq_best[i][j]*yq_best[i][j]);
      }
      for(j=0; j<1*nq; j++)
        fprintf(fpq, "%e %e\n", yq_best[i][j]* dataset[i].con.scale, yq_std[i][j]* dataset[i].con.scale);
      for(j=0; j<dataset[i].nlset; j++)
        for(k=0; k<nq; k++)
          fprintf(fpq, "%e %e\n", yq_best[i][(1+j)*nq+k]* dataset[i].line[j].scale, yq_std[i][(1+j)*nq+k]* dataset[i].line[j].scale);
      fprintf(fpq, "\n");
      
      if(parset.flag_trend > 0)
      {
        printf("Longterm q of dataset %d: Val and Err\n", i);
        for(j=0; j<1*nq; j++)
        {
          printf("%e %e\n", yq_best[i][j]*dataset[i].con.scale, yq_std[i][j]* dataset[i].con.scale);
        }
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nq; k++)
            printf("%e %e\n", yq_best[i][(1+j)*nq+k] * dataset[i].line[j].scale, yq_std[i][(1+j)*nq+k]* dataset[i].line[j].scale);
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
    if(parset.model == pmap)
    {
      free(post_model_trans);
    }
    free(yq);
    free(yq_best);
    free(yq_std);
  }

  return;
}

/*!
 *  this function does postprocess. 
 */
void postprocess_line()
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

/*!
 * matrix operation A^-1 x B is implemented by calling functions
 *    inverse_mat()
 *    multiply_mat_MN()
 */
void reconstruct_line_from_varmodel2(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq;
  int i, j, k, m, info, idx, *ipiv;
  double *PEmat1, *PEmat2, *PEmat3, *PEmat4;
  int nall_data, nqall, ntall, np;
  double *fall_data;
  double sigma, tau, *pm=(double *)model;

  ipiv = workspace_ipiv;
  
  idx = idx_con_pm[nds];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]) * sqrt(tau);

  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  nqall = nq * (1  + dataset[nds].nlset);
  nall_data = alldata[nds].n;
  fall_data = alldata[nds].f;

  Larr = workspace;
  ybuf = Larr + nall_data * nqall;
  y = ybuf + nall_data * nqall;
  Cq = y + nall_data;
  yq = Cq + nqall*nqall;
  yuq = yq + nqall; 
  Larr_rec = yuq + ntall;

  PEmat1 = malloc(ntall * nall_data * sizeof(double));
  PEmat2 = malloc(ntall * ntall * sizeof(double));
  PEmat3 = malloc(ntall * nqall * sizeof(double));
  PEmat4 = malloc(ntall * ntall * sizeof(double));

  set_covar_Pmat_data_line(model, nds);
  
  set_covar_Umat_line(model, nds, nall, tall);
  
  set_covar_Amat_line(model, nds, nall, tall);

  inverse_symat(PCmat, nall_data, &info, ipiv);

  for(i=0;i<dataset[nds].con.n;i++)
  {
    Larr[i*nqall] = 1.0; 
    for(j=1; j<nq; j++)
      Larr[i*nqall+j] = pow(dataset[nds].con.t[i], j);

    for(j=nq; j<nqall; j++)
      Larr[i*nqall + j] = 0.0;
  }
  np = dataset[nds].con.n;
  for(j=0; j<dataset[nds].nlset; j++)
  {
    for(m=0; m<dataset[nds].line[j].n; m++)
    {
      for(i=0; i<nqall; i++)
        Larr[(np+m)*nqall + i ]  = 0.0;
        
      Larr[(np+m)*nqall + nq + j*nq + 0] = 1.0;
      for(i=1; i<nq; i++)
        Larr[(np+m)*nqall + nq + j*nq + i] = pow(dataset[nds].line[j].t[m], i);
    }
    np += dataset[nds].line[j].n;
  }

  /* L^T x C^-1 x L */
  multiply_mat_MN(PCmat, Larr, ybuf, nall_data, nqall, nall_data);
  multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall_data);

  /* L^T x C^-1 x y */
  multiply_matvec(PCmat, fall_data, nall_data, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yuq, nqall, 1, nall_data);

  /* (L^T x C^-1 x L)^-1 x  L^T x C^-1 x y */
  inverse_symat(Cq, nqall, &info, ipiv);
  multiply_mat_MN(Cq, yuq, yq, nqall, 1, nqall);
  memcpy(yqall, yq, nqall*sizeof(double));

  /*  L x q */
  multiply_matvec_MN(Larr, nall_data, nqall, yq, ybuf);
  for(i=0; i<nall_data; i++)
  {
    y[i] = fall_data[i] - ybuf[i];
  }

  /* S x C^-1 x (y - Lq) */
  multiply_matvec(PCmat, y, nall_data, ybuf);
  multiply_matvec_MN(USmat, ntall, nall_data, ybuf, fall);

  for(i=0;i<nall[0];i++)
  {
    Larr_rec[i*nqall + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nqall + j]=pow(tall[i], j);

    for(j=nq; j<nqall; j++)
      Larr_rec[i*nqall + j] = 0.0;
  }

  np = nall[0];
  for(k=0; k<dataset[nds].nlset; k++)
  {
    for(i=0; i<nall[1+k]; i++)
    {
      for(j=0; j<nqall; j++)
      {
        Larr_rec[(np+i)*nqall + j] = 0.0;
      }
      
      Larr_rec[(np+i)*nqall + nq + k*nq + 0] = 1.0;
      for(j=1; j<nq; j++)
        Larr_rec[(np+i)*nqall + nq + k*nq + j] = pow(tall[np+i], j);
    }
    np += nall[1+k];
  }

  multiply_matvec_MN(Larr_rec, ntall, nqall, yq, yuq);

  for(i=0; i<ntall; i++)
    fall[i] += yuq[i];

  /* S x C^-1 x S */
  multiply_mat_MN(USmat, PCmat, PEmat1, ntall, nall_data, nall_data);
  multiply_mat_MN_transposeB(PEmat1, USmat, PEmat2, ntall, ntall, nall_data);

  /* S x C^-1 x L */
  multiply_mat_MN(PEmat1, Larr, PEmat3, ntall, nqall, nall_data);
  /* S x C^-1 x L - L */
  for(i=0; i<ntall*nqall; i++)PEmat3[i] -= Larr_rec[i];
  /* (S x C^-1 x L - L) x Cq */
  multiply_mat_MN(PEmat3, Cq, PEmat1, ntall, nqall, nqall);
  /* (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  multiply_mat_MN_transposeB(PEmat1, PEmat3, PEmat4, ntall, ntall, nqall);

  for(i=0; i<ntall; i++)
  {
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat2[i*ntall+i] + PEmat4[i*ntall + i]);
    //printf("%d %f %f %f\n", i, ASmat[i*ntall + i],  PEmat2[i*ntall+i], PEmat4[i*ntall + i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  return;

}

/*!
 * matrix operation A^-1 x B is implemented by calling function 
 *    multiply_mat_MN_inverseA()
 * 
 */
void reconstruct_line_from_varmodel(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall)
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
  sigma = exp(pm[idx+1]) * sqrt(tau);

  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  nqall = nq * (1  + dataset[nds].nlset);
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

  set_covar_Pmat_data_line(model, nds);
  
  set_covar_Umat_line(model, nds, nall, tall);
  
  set_covar_Amat_line(model, nds, nall, tall);

  for(i=0;i<dataset[nds].con.n;i++)
  {
    Larr[i*nqall] = 1.0; 
    for(j=1; j<nq; j++)
      Larr[i*nqall+j] = pow(dataset[nds].con.t[i], j);

    for(j=nq; j<nqall; j++)
      Larr[i*nqall + j] = 0.0;
  }
  np = dataset[nds].con.n;
  for(j=0; j<dataset[nds].nlset; j++)
  {
    for(m=0; m<dataset[nds].line[j].n; m++)
    {
      for(i=0; i<nqall; i++)
        Larr[(np+m)*nqall + i ]  = 0.0;
        
      Larr[(np+m)*nqall + nq + j*nq + 0] = 1.0;
      for(i=1; i<nq; i++)
        Larr[(np+m)*nqall + nq + j*nq + i] = pow(dataset[nds].line[j].t[m], i);
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

  for(i=0;i<nall[0];i++)
  {
    Larr_rec[i*nqall + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nqall + j]=pow(tall[i], j);

    for(j=nq; j<nqall; j++)
      Larr_rec[i*nqall + j] = 0.0;
  }

  np = nall[0];
  for(k=0; k<dataset[nds].nlset; k++)
  {
    for(i=0; i<nall[1+k]; i++)
    {
      for(j=0; j<nqall; j++)
      {
        Larr_rec[(np+i)*nqall + j] = 0.0;
      }

      Larr_rec[(np+i)*nqall + nq + k*nq + 0] = 1.0;
      for(j=1; j<nq; j++)
        Larr_rec[(np+i)*nqall + nq + k*nq + j] = pow(tall[np+i], j);
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

  /* S x C^-1 x L */
  multiply_mat_MN_transposeA(PEmat2, Larr, PEmat3, ntall, nqall, nall_data);
  /* S x C^-1 x L - L */
  for(i=0; i<ntall*nqall; i++)PEmat3[i] -= Larr_rec[i];
  
  inverse_mat(Cq, nqall, &info, ipiv);

  multiply_mat_MN(PEmat3, Cq, PEmat2, ntall, nqall, nqall);
  /* (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  multiply_mat_MN_transposeB(PEmat2, PEmat3, PEmat4, ntall, ntall, nqall);
  
  for(i=0; i<ntall; i++)
  {
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat1[i*ntall+i] + PEmat4[i*ntall + i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  return;
}

/*!
 * calculate inverse using matrix partition 
 * 
 */
void reconstruct_line_from_varmodel3(const void *model, int nds, int *nall, double *tall, 
                                    double *fall, double *feall, double *yqall)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq, *W, *D, *phi, *fe;
  int i, j, k, m, info, idx, *ipiv;
  double *PEmat1, *PEmat2;
  int nall_data, nqall, ntall, np;
  double *fall_data;
  double sigma, tau, *pm=(double *)model;
  double *work, syserr, lndet, sigma2;
  int *narr, nd;

  ipiv = workspace_ipiv;
  work = workspace_inv;
  
  idx = idx_con_pm[nds];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]) * sqrt(tau);
  /* note here divided by sigma */
  syserr = (exp(pm[idx+0])-1.0)*dataset[nds].con.error_mean/sigma;
  sigma2 = sigma * sigma;

  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  nqall = nq * (1  + dataset[nds].nlset);
  nall_data = alldata[nds].n;
  fall_data = alldata[nds].f;

  Larr = workspace;
  ybuf = Larr + nall_data * nqall;
  y = ybuf + nall_data * nqall;
  Cq = y + nall_data;
  yq = Cq + nqall*nqall;
  yuq = yq + nqall; 
  Larr_rec = yuq + ntall;

  W = Larr_rec + ntall * nqall;
  D = W + nall_data;
  phi = D + nall_data;
  fe = phi + nall_data;
  narr = (int *)(fe + dataset[nds].con.n);

  PEmat1 = malloc(ntall * nall_data * sizeof(double));
  PEmat2 = malloc(ntall * ntall * sizeof(double));

  set_covar_Pmat_data_line(model, nds);
  
  set_covar_Umat_line(model, nds, nall, tall);
  
  set_covar_Amat_line(model, nds, nall, tall);

  for(i=0;i<dataset[nds].con.n;i++)
  {
    Larr[i*nqall] = 1.0; 
    for(j=1; j<nq; j++)
      Larr[i*nqall+j] = pow(dataset[nds].con.t[i], j);

    for(j=nq; j<nqall; j++)
      Larr[i*nqall + j] = 0.0;
  }
  np = dataset[nds].con.n;
  for(j=0; j<dataset[nds].nlset; j++)
  {
    for(m=0; m<dataset[nds].line[j].n; m++)
    {
      for(i=0; i<nqall; i++)
        Larr[(np+m)*nqall + i ]  = 0.0;
        
      Larr[(np+m)*nqall + nq + j*nq + 0] = 1.0;
      for(i=1; i<nq; i++)
        Larr[(np+m)*nqall + nq + j*nq + i] = pow(dataset[nds].line[j].t[m], i);
    }
    np += dataset[nds].line[j].n;
  }
  
  nd = dataset[nds].nlset;
  narr[0] = dataset[nds].con.n;
  for(i=1; i<=nd; i++)
    narr[i] = dataset[nds].line[i-1].n;
  
  for(i=0; i<dataset[nds].con.n; i++)
    fe[i] = dataset[nds].con.fe[i]/sigma;

  /* note here set sigma = 1 */
  inverse_semiseparable_iter(dataset[nds].con.t, dataset[nds].con.n, 1.0, 1.0/tau, 
                          fe, syserr, W, D, phi, PCmat, nall_data, narr, nd,
                          &lndet, work, ipiv);

  /* L^T x C^-1 x L */
  multiply_mat_MN(PCmat, Larr, ybuf, nall_data, nqall, nall_data);
  multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall_data);

  /* L^T x C^-1 x y, note ybuf = C^-1 x L */
  multiply_matvec_MN_transposeA(ybuf, nall_data, nqall, fall_data, yuq);

  /* (L^T x C^-1 x L)^-1 x  L^T x C^-1 x y */
  inverse_symat(Cq, nqall, &info, ipiv);
  multiply_mat_MN(Cq, yuq, yq, nqall, 1, nqall);
  for(i=0; i<nqall; i++)
  {
    yqall[i] = yq[i];
  }

  /*  y = -1.0 * L x q + 1.0 * y */
  memcpy(y, fall_data, nall_data*sizeof(double));
  multiply_matvec_MN_alpha_beta(Larr, nall_data, nqall, yq, y, -1.0, 1.0);

  /* S x C^-1 x (y - Lq) */  
  /* PEmat1 = S x C^-1, dimension: n*nd */
  multiply_mat_MN(USmat, PCmat, PEmat1, ntall, nall_data, nall_data);
  multiply_matvec_MN(PEmat1, ntall, nall_data, y, fall);

  for(i=0;i<nall[0];i++)
  {
    Larr_rec[i*nqall + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nqall + j]=pow(tall[i], j);

    for(j=nq; j<nqall; j++)
      Larr_rec[i*nqall + j] = 0.0;
  }

  np = nall[0];
  for(k=0; k<dataset[nds].nlset; k++)
  {
    for(i=0; i<nall[1+k]; i++)
    {
      for(j=0; j<nqall; j++)
      {
        Larr_rec[(np+i)*nqall + j] = 0.0;
      }
      
      Larr_rec[(np+i)*nqall + nq + k*nq + 0] = 1.0;
      for(j=1; j<nq; j++)
        Larr_rec[(np+i)*nqall + nq + k*nq + j] = pow(tall[np+i], j);
    }
    np += nall[1+k];
  }
  
  /* fall + L x q */
  multiply_matvec_MN_alpha_beta(Larr_rec, ntall, nqall, yq, fall, 1.0, 1.0);

  /* S x C^-1 x S */
  multiply_mat_MN_transposeB(PEmat1, USmat, PEmat2, ntall, ntall, nall_data);
  /* S x C^-1 x L - L */
  multiply_mat_MN_alpha_beta(PEmat1, Larr, Larr_rec, ntall, nqall, nall_data, 1.0, -1.0);
  /* (S x C^-1 x L - L) x Cq */
  multiply_mat_MN(Larr_rec, Cq, PEmat1, ntall, nqall, nqall);
  /* -S x C^-1 x S + (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  multiply_mat_MN_transposeB_alpha_beta(PEmat1, Larr_rec, PEmat2, ntall, ntall, nqall, 1.0, -1.0);

  for(i=0; i<ntall; i++)
  {
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] + PEmat2[i*ntall+i]);
    //printf("%d %f %f %f\n", i, ASmat[i*ntall + i],  PEmat2[i*ntall+i], PEmat4[i*ntall + i]);
  }

  free(PEmat1);
  free(PEmat2);
  return;
}

/*!
 * matrix operation A^-1 x B is implemented by calling function 
 *    multiply_mat_MN_inverseA()
 * 
 */
double prob_line_variability(const void *model)
{
  double prob = 0.0, prob1;
  int i, j, k, m, np, info, sign, *ipiv;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq, *ICq, *yave;
  double *fall;
  int nall, nqall, idx;
  double sigma, tau, *pm = (double *)model;

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max;
  yq = y + nall_max;
  yave = yq + (1+nlset_max)*nq;
  Cq = yave + nall_max;
  ICq = Cq + (1+nlset_max)*nq * (1+nlset_max)*nq;

  ipiv = workspace_ipiv;

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    idx = idx_con_pm[k];
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]) * sqrt(tau);

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (1+dataset[k].nlset);

    for(i=0;i<dataset[k].con.n;i++)
    {
      Larr[i*nqall] = 1.0; 
      for(j=1; j<nq; j++)
        Larr[i*nqall + j] = pow(dataset[k].con.t[i], j);

      for(j=nq; j<nqall; j++)
        Larr[i*nqall + j] = 0.0;
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
     
    set_covar_Pmat_data_line_array(model, k);
    
    memcpy(Tmat1, PCmat, nall*nall*sizeof(double));
    memcpy(Tmat2, Larr, nall*nqall*sizeof(double));

    multiply_mat_MN_inverseA(Tmat1, Tmat2, nall, nqall, ipiv); // Tmat2 = C^-1 * L;  NxNq

    multiply_mat_MN_transposeA(Larr, Tmat2, ICq, nqall, nqall, nall); // ICq = L^T*C^-1*L; NqxNq
    multiply_mat_MN_transposeA(Tmat2, fall, yq, nqall, 1, nall); // yq = L^T*C^-1*y;  Nqx1
    memcpy(Tmat1, ICq, nqall*nqall*sizeof(double));
    multiply_mat_MN_inverseA(Tmat1, yq, nqall, 1, ipiv); // yq = (L^T*C^-1*L)^-1 * L^T*C^-1*y; Nqx1

    multiply_mat_MN(Larr, yq, yave, nall, 1, nqall); // yave = L * q; Nx1
    for(i=0; i<nall; i++)y[i] = fall[i] - yave[i];
    memcpy(Tmat1, PCmat, nall*nall*sizeof(double));
    memcpy(ybuf, y, nall*sizeof(double));
    multiply_mat_MN_inverseA(Tmat1, ybuf, nall, 1, ipiv); // ybuf = C^-1 * y; Nx1

    prob1 = -0.5*cblas_ddot(nall, y, 1, ybuf, 1)/(sigma*sigma); // y^T * C^-1 * y

    if(prob1 > 0.0 )  // check if prob is positive
    { 
      prob = -DBL_MAX;
      printf("prob >0!\n");
      return prob;
    }
    lndet = lndet_mat3(PCmat, nall, &info, &sign, ipiv) + 2.0*nall*log(sigma);
    if(info!=0|| sign==-1)
    {
      prob = -DBL_MAX;
      printf("lndet_C %f %d!\n", lndet, sign);
      return prob;
    }
    lndet_ICq = lndet_mat3(ICq, nqall, &info, &sign, ipiv) - 2.0*nqall*log(sigma);
    if(info!=0 || sign==-1 )
    {
      prob = -DBL_MAX;
      printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
      return prob;
    }
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;
  }

  return prob;
}

/*!
 * matrix operation A^-1 x B is implemented by calling functions
 *    inverse_mat()
 *    multiply_mat_MN()
 */
double prob_line_variability2(const void *model)
{
  double prob = 0.0, prob1, sigma, tau;
  int i, j, k, m, np, info, sign, *ipiv;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq, *ICq;
  double *fall;
  int nall, nqall, idx;
  double *pm = (double *)model;

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max * ((1+nlset_max)*nq);
  yq = y + nall_max;
  Cq = yq + (1+nlset_max)*nq;
  ICq = Cq + (1+nlset_max)*nq * (1+nlset_max)*nq;

  ipiv = workspace_ipiv;

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    idx = idx_con_pm[k];
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]) * sqrt(tau);

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (1+dataset[k].nlset);

    for(i=0;i<dataset[k].con.n;i++)
    {
      Larr[i*nqall] = 1.0; 
      for(j=1; j<nq; j++)
        Larr[i*nqall + j] = pow(dataset[k].con.t[i], j);

      for(j=nq; j<nqall; j++)
        Larr[i*nqall + j] = 0.0;
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
     
    set_covar_Pmat_data_line_array(model, k);
    memcpy(IPCmat, PCmat, nall*nall*sizeof(double));

    inverse_mat(IPCmat, nall, &info, ipiv); /* calculate C^-1 */
    
    /* calculate L^T*C^-1*L */
    multiply_mat_MN(IPCmat, Larr, ybuf, nall, nqall, nall);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall);
    memcpy(ICq, Cq, nqall*nqall*sizeof(double));

    /* calculate L^T*C^-1*y */
    multiply_matvec(IPCmat, fall, nall, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nqall, 1, nall);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_mat(Cq, nqall, &info, ipiv);
    multiply_mat_MN(Cq, yq, ybuf, nqall, 1, nqall);
  
    multiply_matvec_MN(Larr, nall, nqall, ybuf, y);
    for(i=0; i<nall; i++)
    {
      y[i] = fall[i] - y[i];
    }

    /* y^T x C^-1 x y */
    multiply_matvec(IPCmat, y, nall, ybuf);
    prob1 = -0.5 * cblas_ddot(nall, y, 1, ybuf, 1)/(sigma*sigma);
    
    if(prob1 > 0.0 )  // check if prob is positive
    { 
      prob = -DBL_MAX;
      printf("prob >0!\n");
      return prob;
    }
    lndet = lndet_mat3(PCmat, nall, &info, &sign, ipiv) + 2.0*nall*log(sigma);
    if(info!=0|| sign==-1)
    {
      prob = -DBL_MAX;
      printf("lndet_C %f %d!\n", lndet, sign);
      return prob;
    }
    lndet_ICq = lndet_mat3(ICq, nqall, &info, &sign, ipiv) - 2.0*nqall*log(sigma);
    if(info!=0 || sign==-1 )
    {
      prob = -DBL_MAX;
      printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
      return prob;
    }
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;

  }
  return prob;
}

/*!
 * calulate inverse and lndet simultaneously.
 */
double prob_line_variability3(const void *model)
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
    sigma = exp(pm[idx+1]) * sqrt(tau);

    nall = alldata[k].n;
    fall = alldata[k].f;
    nqall = nq * (1+dataset[k].nlset);

    for(i=0;i<dataset[k].con.n;i++)
    {
      Larr[i*nqall] = 1.0; 
      for(j=1; j<nq; j++)
        Larr[i*nqall + j] = pow(dataset[k].con.t[i], j);

      for(j=nq; j<nqall; j++)
        Larr[i*nqall + j] = 0.0;
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
     
    set_covar_Pmat_data_line_array(model, k);

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
 * using semiseparable matrix and partitioning inverse
 * 
 */
double prob_line_variability4(const void *model)
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

    set_covar_Pmat_data_line_array(model, k);

    for(i=0; i<dataset[k].con.n; i++)
      fe[i] = dataset[k].con.fe[i]/sigma;
      
    inverse_semiseparable_iter(dataset[k].con.t, dataset[k].con.n, 1.0, 1.0/tau, 
                          fe, syserr, W, D, phi, PCmat, nall, narr, nd,
                          &lndet, work, ipiv);

    /* C^-1 and lndet(C) */
    //inverse_symat_partition_iter(PCmat, nall, narr, nd, &lndet, work, ipiv); 
    lndet += 2.0*nall*log(sigma);

    /* calculate L^T*C^-1*L, ybuf = C^-1xL */
    multiply_mat_MN(PCmat, Larr, ybuf, nall, nqall, nall);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall);

    /* calculate L^T*C^-1*y */
    // multiply_matvec(PCmat, fall, nall, ybuf);
    // multiply_mat_MN_transposeA(Larr, ybuf, yq, nqall, 1, nall);

    /* L^T*C^-1*y = (C^-1xL)^T x y = ybuf^T x y */
    multiply_matvec_MN_transposeA(ybuf, nall, nqall, fall, yq);

    /* calculate q = (L^T*C^-1*L)^-1 * L^T*C^-1*y */
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

/*!
 * check whether q falls within the flux limit
 */
double check_flux(double *yq, int iset)
{
  int j;
  double diff, prob=0.0;

  /*check whether q is beyond flux min and max */
  diff = fabs(yq[0*nq+0] - flux_minmax[iset][0*2+0])/flux_minmax[iset][0*2+1];
  if(diff > 1.0)
    prob -= 100.0*pow(diff, 2);
  for(j=0; j<dataset[iset].nlset; j++)
  {
    diff = fabs(yq[(1+j)*nq+0] - flux_minmax[iset][(1+j)*2+0])/flux_minmax[iset][(1+j)*2+1];
    if(diff > 1.0)
      prob -= 100.0*pow(diff, 2);
  }

  return prob;
}

int mc_line_init()
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
    lag_prior = (double *)malloc(parset.num_gaussian_upper * 2 * sizeof(double));
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

int mc_line_end()
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
 * this function sets the covariance matrix at data time points 
 *
 * k is the index of dataset
 */
void set_covar_Pmat_data_line(const void *model, int k)
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
  sigma = exp(pm[idx+1]) * sqrt(tau);
  syserr = (exp(pm[idx+0])-1.0)*dataset[k].con.error_mean;
  

  /* continuum - continuum  */
  for(i=0; i<dataset[k].con.n; i++)
  { 
    t1 = dataset[k].con.t[i];

    for(j=0; j<i; j++)
    {
      t2 = dataset[k].con.t[j];
      PCmat[i*nall+j] = PCmat[j*nall+i] = 
      PSmat[i*nall+j] = PSmat[j*nall+i] = exp (- fabs(t1-t2) / tau );
    }
    PSmat[i*nall+i] = 1.0;
    error = dataset[k].con.fe[i]*dataset[k].con.fe[i] + syserr*syserr;
    PCmat[i*nall+i] = PSmat[i*nall+i] + error/sigma/sigma;
  }

  /* continuum - line */
  for(i=0; i<dataset[k].con.n; i++)
  { 
    t1 = dataset[k].con.t[i];
    np = dataset[k].con.n;
    for(j=0; j<dataset[k].nlset; j++)
    {
      for(m=0; m < dataset[k].line[j].n; m++)
      {
        t2 = dataset[k].line[j].t[m];
        PCmat[i*nall + np + m] = PCmat[(np+m)*nall + i] =
        PSmat[i*nall + np + m] = PSmat[(np+m)*nall + i] = Slc(t1, t2, model, k, j);
      }
      np += dataset[k].line[j].n;
    }
  }

  /* line - line */
  np = dataset[k].con.n;
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
 * this function sets the covariance matrix at data time points for an array of times 
 *
 * k is the index of dataset
 */
void set_covar_Pmat_data_line_array(const void *model, int k)
{
  double t1, t2;
  double sigma, tau, syserr, syserr_line, error;
  int i, j, m, l, nall, np, npline, idx, ncon, nline, nline2;
  double *pm = (double *)model, *tcon, *tline, *tline2;
  
  /* total number of point of this dataset */
  nall = alldata[k].n;

  /* set variability parameters */
  idx = idx_con_pm[k];
  tau = exp(pm[idx+2]);
  sigma = exp(pm[idx+1]) * sqrt(tau);
  syserr = (exp(pm[idx+0])-1.0)*dataset[k].con.error_mean;
  

  /* continuum - continuum  */
  for(i=0; i<dataset[k].con.n; i++)
  { 
    t1 = dataset[k].con.t[i];

    for(j=0; j<i; j++)
    {
      t2 = dataset[k].con.t[j];
      PCmat[i*nall+j] = PCmat[j*nall+i] = 
      PSmat[i*nall+j] = PSmat[j*nall+i] = exp (- fabs(t1-t2) / tau );
    }
    PSmat[i*nall+i] = 1.0;
    error = dataset[k].con.fe[i]*dataset[k].con.fe[i] + syserr*syserr;
    PCmat[i*nall+i] = PSmat[i*nall+i] + error/sigma/sigma;
  }

  /* continuum - line */
  tcon = dataset[k].con.t;
  ncon = dataset[k].con.n;
  np = ncon;
  for(j=0; j<dataset[k].nlset; j++)
  {
    tline = dataset[k].line[j].t;
    nline = dataset[k].line[j].n;
    Slc_array(tcon, ncon, tline, nline, model, k, j, Smat_lc);
    for(i=0; i<ncon; i++)
    {
      for(m=0; m < nline; m++)
      {
        PCmat[i*nall + np + m] = PCmat[(np+m)*nall + i] =
        PSmat[i*nall + np + m] = PSmat[(np+m)*nall + i] = Smat_lc[i*nline + m];
      }
    }
    np += nline;
  }

  /* line - line */
  np = ncon;
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

/*!
 * this function sets the covariance matrix at time of data points and reconstruction points
 */
void set_covar_Umat_line(const void *model, int nds, int *nall, double *tall)
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
    t1 = tall[i];
   
    /* continuum self */
    for(j=0; j<dataset[nds].con.n; j++)
    {
      t2 = dataset[nds].con.t[j];
      USmat[i * nall_data + j] =  exp( - fabs(t1-t2)/taud);
    }

    /* lines */
    np = dataset[nds].con.n;
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
  
  /* line - continuum/line */
  npline = nall[0];
  for(k = 0; k<dataset[nds].nlset; k++)
  {
    
    for(i=0; i<nall[1+k]; i++)
    {
      t1 = tall[npline + i];

      /* continum */
      for(j=0; j<dataset[nds].con.n; j++)
      {
        t2 = dataset[nds].con.t[j];
        USmat[ (i+npline) * nall_data + j ] = Slc(t2, t1, model, nds, k);
      }

      /* line different */
      np = dataset[nds].con.n;
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
 */
void set_covar_Amat_line(const void *model, int nds, int *nall, double *tall)
{
  double t1, t2, taud;
  int i, j, k, m, ntall, np, npline, idx;
  double *pm=(double *)model;
  
  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  /* continuum - continuum/lines*/
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

/*===============================================================================
 * Gaussian transfer function 
 *===============================================================================
 */
/*
 * auto-covariance of line
 *
 * nds: index of dataset
 * nls: index of line
 *
 */
double Sll_gauss(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot, A;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      DT = Dt - (tau1 - tau2);
    
     /* 
      * for very large x 
      * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
      * 
      * (see https://en.wikipedia.org/wiki/Error_function)
      *
      * with this equation, the factor exp(wg*wg/taud/taud) will be canceled out.
      */
     /*
      {
        St = exp(wg*wg/taud/taud) * ( exp(-DT/taud) * erfc( -DT/2.0/wg + wg/taud )
                                 +exp( DT/taud) * erfc(  DT/2.0/wg + wg/taud ) );
      }
     */

      A = sqrt(wg1*wg1 + wg2*wg2);

      St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud )
          +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud ) ;

      St *= 1.0/2.0 * fg1*fg2;

      Sttot += St;
    }
  }
  
  return Sttot;
}

double Sll_gauss_linear(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot, A;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      DT = Dt - (tau1 - tau2);
    
     /* 
      * for very large x 
      * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
      * 
      * (see https://en.wikipedia.org/wiki/Error_function)
      *
      * with this equation, the factor exp(wg*wg/taud/taud) will be canceled out.
      */
     /*
      {
        St = exp(wg*wg/taud/taud) * ( exp(-DT/taud) * erfc( -DT/2.0/wg + wg/taud )
                                 +exp( DT/taud) * erfc(  DT/2.0/wg + wg/taud ) );
      }
     */

      A = sqrt(wg1*wg1 + wg2*wg2);

      St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud )
          +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud ) ;

      St *= 1.0/2.0 * fg1*fg2;

      Sttot += St;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of line for an array of times
 *
 * nds: index of dataset
 * nls: index of line
 *
 */

void Sll_array_gauss(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St, A, At, At2;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, fg12;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      A = sqrt(wg1*wg1 + wg2*wg2);
      At = A/taud;
      At2 = At*At/2.0;

      fg12 = fg1*fg2/2.0;

      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau1 - tau2);

          St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - At) / sqrt(2.0) ) + At2 )
              +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + At) / sqrt(2.0) ) + At2 ) ;

          St *= fg12;

          Smat[i*nline + j] += St;
        }
      }
    }
  }

  return;
}

void Sll_array_gauss_linear(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St, A, At, At2;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, fg12;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      A = sqrt(wg1*wg1 + wg2*wg2);
      At = A/taud;
      At2 = At*At/2.0;

      fg12 = fg1*fg2/2.0;

      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau1 - tau2);

          St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - At) / sqrt(2.0) ) + At2 )
              +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + At) / sqrt(2.0) ) + At2 ) ;

          St *= fg12;

          Smat[i*nline + j] += St;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different line
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_gauss(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau1, tau2, wg1, wg2, taud;
  double Dt, DT, St, Sttot, A;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

  
      DT = Dt - (tau1-tau2);
  
      A = sqrt(wg1*wg1 + wg2*wg2);

      St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud )
          +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud ) ;

      St *= 1.0/2.0 * fg1*fg2;

      Sttot += St;
    }
  }

  return Sttot;
}

double Sll2_gauss_linear(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau1, tau2, wg1, wg2, taud;
  double Dt, DT, St, Sttot, A;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

  
      DT = Dt - (tau1-tau2);
  
      A = sqrt(wg1*wg1 + wg2*wg2);

      St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud )
          +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + A/taud) / sqrt(2.0) ) + A*A/2.0/taud/taud ) ;

      St *= 1.0/2.0 * fg1*fg2;

      Sttot += St;
    }
  }

  return Sttot;
}


/*
 * covariance between different line for an array of times
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */

void Sll2_array_gauss(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau1, tau2, wg1, wg2, taud;
  double Dt, DT, St, A, At, At2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      A = sqrt(wg1*wg1 + wg2*wg2);
      At = A/taud;
      At2 = At*At/2.0;

      fg12 = fg1*fg2/2.0;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau1-tau2);

          St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - At) / sqrt(2.0) ) + At2 )
              +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + At) / sqrt(2.0) ) + At2 ) ;

          St *= fg12;

          Smat[i*nline2 + j] += St;
        }
      }
    }
  }

  return;
}

void Sll2_array_gauss_linear(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau1, tau2, wg1, wg2, taud;
  double Dt, DT, St, A, At, At2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      A = sqrt(wg1*wg1 + wg2*wg2);
      At = A/taud;
      At2 = At*At/2.0;

      fg12 = fg1*fg2/2.0;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau1-tau2);

          St = exp( -DT/taud + gsl_sf_log_erfc( -(DT/A - At) / sqrt(2.0) ) + At2 )
              +exp(  DT/taud + gsl_sf_log_erfc(  (DT/A + At) / sqrt(2.0) ) + At2 ) ;

          St *= fg12;

          Smat[i*nline2 + j] += St;
        }
      }
    }
  }

  return;
}

/*
 * covariance between continuum and line 
 *
 * nds: dataset index
 * nls: line set index
 */
double Slc_gauss(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, Sttot;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);

    DT = Dt - tau0;

  /* 
   * for very large x 
   * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
   * 
   * (see https://en.wikipedia.org/wiki/Error_function)
   *
   * with this equation, the factor exp(wg*wg/taud/taud/2) will be canceled out.
   */
  /*
  {
    St = exp(wg*wg/2.0/taud/taud) * ( exp(-DT/taud) * erfc( -(DT/wg - wg/taud)/sqrt(2.0) ) 
                                     +exp( DT/taud) * erfc(  (DT/wg + wg/taud)/sqrt(2.0) ));
  }
  */

    St = exp(-DT/taud + gsl_sf_log_erfc( -(DT/wg - wg/taud)/sqrt(2.0) ) +  wg*wg/2.0/taud/taud )
        +exp( DT/taud + gsl_sf_log_erfc(  (DT/wg + wg/taud)/sqrt(2.0) ) +  wg*wg/2.0/taud/taud );

    St *= 1.0/2.0 * fg;

    Sttot += St;
  }

  return Sttot;
}

double Slc_gauss_linear(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, Sttot;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);

    DT = Dt - tau0;

  /* 
   * for very large x 
   * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
   * 
   * (see https://en.wikipedia.org/wiki/Error_function)
   *
   * with this equation, the factor exp(wg*wg/taud/taud/2) will be canceled out.
   */
  /*
  {
    St = exp(wg*wg/2.0/taud/taud) * ( exp(-DT/taud) * erfc( -(DT/wg - wg/taud)/sqrt(2.0) ) 
                                     +exp( DT/taud) * erfc(  (DT/wg + wg/taud)/sqrt(2.0) ));
  }
  */

    St = exp(-DT/taud + gsl_sf_log_erfc( -(DT/wg - wg/taud)/sqrt(2.0) ) +  wg*wg/2.0/taud/taud )
        +exp( DT/taud + gsl_sf_log_erfc(  (DT/wg + wg/taud)/sqrt(2.0) ) +  wg*wg/2.0/taud/taud );

    St *= 1.0/2.0 * fg;

    Sttot += St;
  }

  return Sttot;
}

/*
 * covariance between continuum and line for an array of times
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_gauss(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, wt, wt2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    wt = wg/taud;
    wt2 = wt*wt/2.0;
    
    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        St = exp(-DT/taud + gsl_sf_log_erfc( -(DT/wg - wt)/sqrt(2.0) ) +  wt2 )
            +exp( DT/taud + gsl_sf_log_erfc(  (DT/wg + wt)/sqrt(2.0) ) +  wt2 );

        St *= 1.0/2.0 * fg;

        Smat[i*nline + j] += St;
      }
    }
  }

  return;
}

void Slc_array_gauss_linear(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, wt, wt2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    wt = wg/taud;
    wt2 = wt*wt/2.0;
    
    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        St = exp(-DT/taud + gsl_sf_log_erfc( -(DT/wg - wt)/sqrt(2.0) ) +  wt2 )
            +exp( DT/taud + gsl_sf_log_erfc(  (DT/wg + wt)/sqrt(2.0) ) +  wt2 );

        St *= 1.0/2.0 * fg;

        Smat[i*nline + j] += St;
      }
    }
  }

  return;
}

/*===============================================================================
 * tophat transfer function 
 *===============================================================================
 */
/*
 * covariance between continuum and line
 *
 * nds: dataset index
 * nls: line set index
 */
double Slc_tophat(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, Sttot;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    fg *= taud/2.0/wg;

    DT = Dt - tau0;

    if(DT < -wg)
    {
      St = exp( (DT + wg) / taud ) - exp( (DT - wg)/taud);
    }
    else if(DT < wg)
    {
      St = 2.0 - exp(- (DT + wg) / taud) - exp( (DT - wg)/taud);
    }
    else 
    {
      St = exp( - (DT - wg)/taud) - exp ( - (DT + wg)/taud);
    }

    Sttot += St * fg;
  }

  return Sttot;
}

double Slc_tophat_linear(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St, Sttot;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    fg *= taud/2.0/wg;

    DT = Dt - tau0;

    if(DT < -wg)
    {
      St = exp( (DT + wg) / taud ) - exp( (DT - wg)/taud);
    }
    else if(DT < wg)
    {
      St = 2.0 - exp(- (DT + wg) / taud) - exp( (DT - wg)/taud);
    }
    else 
    {
      St = exp( - (DT - wg)/taud) - exp ( - (DT + wg)/taud);
    }

    Sttot += St * fg;
  }

  return Sttot;
}

/*
 * covariance between continuum and line for an array of times
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_tophat(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    fg *= taud/2.0/wg;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT < -wg)
        {
          St = exp( (DT + wg) / taud ) - exp( (DT - wg)/taud);
        }
        else if(DT < wg)
        {
          St = 2.0 - exp(- (DT + wg) / taud) - exp( (DT - wg)/taud);
        }
        else 
        {
          St = exp( - (DT - wg)/taud) - exp ( - (DT + wg)/taud);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

void Slc_array_tophat_linear(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, wg, tau0, St;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    wg = exp(pm[idxk + 2]);
    
    fg *= taud/2.0/wg;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT < -wg)
        {
          St = exp( (DT + wg) / taud ) - exp( (DT - wg)/taud);
        }
        else if(DT < wg)
        {
          St = 2.0 - exp(- (DT + wg) / taud) - exp( (DT - wg)/taud);
        }
        else 
        {
          St = exp( - (DT - wg)/taud) - exp ( - (DT + wg)/taud);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

/*
 * auto-covariance between a line for an array of times
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
double Sll_tophat(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, wh, wl, fg12;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }

      DT = Dt - (tau1 - tau2);
    
      if(DT < -wh-wl)
      {
        St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
      }
      else if(DT < -wh+wl)
      {        
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
            +2.0*(DT+wh+wl)/taud;
      }
      else if(DT <  wh-wl)
      {
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
            +4.0*wl/taud;
      }
      else if(DT <  wh+wl)
      {
        St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud)
            -2.0*(DT-wh-wl)/taud;
      }
      else 
      {
        St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

double Sll_tophat_linear(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, wh, wl, fg12;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }

      DT = Dt - (tau1 - tau2);
    
      if(DT < -wh-wl)
      {
        St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
      }
      else if(DT < -wh+wl)
      {        
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
            +2.0*(DT+wh+wl)/taud;
      }
      else if(DT <  wh-wl)
      {
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
            +4.0*wl/taud;
      }
      else if(DT <  wh+wl)
      {
        St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud)
            -2.0*(DT-wh-wl)/taud;
      }
      else 
      {
        St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of a line for an array of times
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
void Sll_array_tophat(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, fg12, wh, wl;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 *taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau1 - tau2);

          if(DT < -wh-wl)
          {
            St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
          }
          else if(DT < -wh+wl)
          {            
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
                +2.0*(DT+wh+wl)/taud;
          }
          else if(DT <  wh-wl)
          {
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                +4.0*wl/taud;
          }
          else if(DT <  wh+wl)
          {   
            St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud) 
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                -2.0*(DT-wh-wl)/taud;
          }
          else 
          {
            St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

void Sll_array_tophat_linear(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau1, wg1, fg2, tau2, wg2, fg12, wh, wl;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 *taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau1 - tau2);

          if(DT < -wh-wl)
          {
            St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
          }
          else if(DT < -wh+wl)
          {            
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
                +2.0*(DT+wh+wl)/taud;
          }
          else if(DT <  wh-wl)
          {
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                +4.0*wl/taud;
          }
          else if(DT <  wh+wl)
          {   
            St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud) 
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                -2.0*(DT-wh-wl)/taud;
          }
          else 
          {
            St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different lines
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_tophat(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau1, tau2, wg1, wg2, taud, wh, wl, fg12;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }
  
      DT = Dt - (tau1-tau2);
  
      if(DT < -wh-wl)
      {
        St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
      }
      else if(DT < -wh+wl)
      {        
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
            +2.0*(DT+wh+wl)/taud;
      }
      else if(DT <  wh-wl)
      {
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
            +4.0*wl/taud;
      }
      else if(DT <  wh+wl)
      {
        St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud)
            -2.0*(DT-wh-wl)/taud;
      }
      else 
      {
        St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

double Sll2_tophat_linear(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau1, tau2, wg1, wg2, taud, wh, wl, fg12;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }
  
      DT = Dt - (tau1-tau2);
  
      if(DT < -wh-wl)
      {
        St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
      }
      else if(DT < -wh+wl)
      {        
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
            +2.0*(DT+wh+wl)/taud;
      }
      else if(DT <  wh-wl)
      {
        St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
            +4.0*wl/taud;
      }
      else if(DT <  wh+wl)
      {
        St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud)
            -2.0*(DT-wh-wl)/taud;
      }
      else 
      {
        St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
            -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

/*
 * covariance between different lines for an array of times
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
void Sll2_array_tophat(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau1, tau2, wg1, wg2, taud, wh, wl;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau1-tau2);

          if(DT < -wh-wl)
          {
            St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
          }
          else if(DT < -wh+wl)
          {            
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
                +2.0*(DT+wh+wl)/taud;
          }
          else if(DT <  wh-wl)
          {
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                +4.0*wl/taud;
          }
          else if(DT <  wh+wl)
          {   
            St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud) 
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                -2.0*(DT-wh-wl)/taud;
          }
          else 
          {
            St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
          }

          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

void Sll2_array_tophat_linear(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau1, tau2, wg1, wg2, taud, wh, wl;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau1 =    pm[idxk1 + 1] ;
    wg1 = exp(pm[idxk1 + 2]);

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau2 =    pm[idxk2 + 1] ;
      wg2 = exp(pm[idxk2 + 2]);

      fg12 = fg1/wg1*fg2/wg2/4.0 * taud*taud;

      if(wg1>=wg2)
      {
        wh = wg1;
        wl = wg2;
      }
      else
      {
        wh = wg2;
        wl = wg1;
      }

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau1-tau2);

          if(DT < -wh-wl)
          {
            St = exp((DT+wh+wl)/taud) + exp((DT-wh-wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud);
          }
          else if(DT < -wh+wl)
          {            
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT+wh-wl)/taud) - exp((DT-wh+wl)/taud) 
                +2.0*(DT+wh+wl)/taud;
          }
          else if(DT <  wh-wl)
          {
            St = exp((DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp((DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                +4.0*wl/taud;
          }
          else if(DT <  wh+wl)
          {   
            St = exp(-(DT+wh+wl)/taud) + exp((DT-wh-wl)/taud) 
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud) 
                -2.0*(DT-wh-wl)/taud;
          }
          else 
          {
            St = exp(-(DT-wh-wl)/taud) + exp(-(DT+wh+wl)/taud)
                -exp(-(DT-wh+wl)/taud) - exp(-(DT+wh-wl)/taud);
          }

          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

/*===============================================================================
 * Gamma transfer function with k=2
 * Gamma(x) = A/T^2 * (x-t0) * exp[-(x-t0)/T]
 * parameters: A, t0, T
 * constraints: x>t0, T>0
 *===============================================================================
 */
/*
 * covaraince between continuum and line
 *
 */
double Slc_gamma(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, Sttot, p1, p2;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);

    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    DT = Dt - tau0;

    if(DT<=0)
    {
      St = exp(DT/taud)/p1/p1;
    }
    else 
    {
      if(fabs(p2) < EPS)
        St = exp(-DT/taud) * ( (DT/taud)*(DT/taud)/2 + (p1*DT/taud+1)/p1/p1 );
      else
        St = exp(-DT/taud)/p2/p2 + exp(-DT/tau1)*(-(p2*DT/tau1+1)/p2/p2 + (p1*DT/tau1+1)/p1/p1);
    }
    
    Sttot += St * fg;
  }

  return Sttot;
}

/*
 * covaraince between continuum and line with linear amplitude
 *
 */
double Slc_gamma_linear(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, Sttot, p1, p2;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);

    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    DT = Dt - tau0;

    if(DT<=0)
    {
      St = exp(DT/taud)/p1/p1;
    }
    else 
    {
      if(fabs(p2) < EPS)
        St = exp(-DT/taud) * ( (DT/taud)*(DT/taud)/2 + (p1*DT/taud+1)/p1/p1 );
      else
        St = exp(-DT/taud)/p2/p2 + exp(-DT/tau1)*(-(p2*DT/tau1+1)/p2/p2 + (p1*DT/tau1+1)/p1/p1);
    }
    
    Sttot += St * fg;
  }

  return Sttot;
}

/*
 * covariance between continuum and line for an array of time
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_gamma(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, p1, p2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);
    
    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT<=0)
        {
          St = exp(DT/taud)/p1/p1;
        }
        else 
        {
          if(fabs(p2) < EPS)
            St = exp(-DT/taud) * ( (DT/taud)*(DT/taud)/2 + (p1*DT/taud+1)/p1/p1 );
          else
            St = exp(-DT/taud)/p2/p2 + exp(-DT/tau1)*(-(p2*DT/tau1+1)/p2/p2 + (p1*DT/tau1+1)/p1/p1);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

/*
 * covariance between continuum and line for an array of time with linear amplitude
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_gamma_linear(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, p1, p2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);
    
    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT<=0)
        {
          St = exp(DT/taud)/p1/p1;
        }
        else 
        {
          if(fabs(p2) < EPS)
            St = exp(-DT/taud) * ( (DT/taud)*(DT/taud)/2 + (p1*DT/taud+1)/p1/p1 );
          else
            St = exp(-DT/taud)/p2/p2 + exp(-DT/tau1)*(-(p2*DT/tau1+1)/p2/p2 + (p1*DT/tau1+1)/p1/p1);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

/*
 * auto-covariance of a line for an array of time
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
double Sll_gamma(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01 - tau02);
    
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
        {
          tmp = p2*DT/tau1;
          St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
        }
        else
          St = exp( DT/taud)/(p2*p2 * p4*p4) 
            + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
      }
      else 
      {
        if(fabs(p1) < EPS)
        {
          tmp = p3*DT/tau2;
          St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
        }
        else 
          St = exp(-DT/taud)/(p1*p1 * p3*p3) 
            + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of a line for an array of time with linear amplitude
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
double Sll_gamma_linear(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =      pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =      pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01 - tau02);
    
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
        {
          tmp = p2*DT/tau1;
          St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
        }
        else
          St = exp( DT/taud)/(p2*p2 * p4*p4) 
            + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
      }
      else 
      {
        if(fabs(p1) < EPS)
        {
          tmp = p3*DT/tau2;
          St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
        }
        else 
          St = exp(-DT/taud)/(p1*p1 * p3*p3) 
            + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of a line for an array of times
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
void Sll_array_gamma(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau01 - tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
            {
              tmp = p2*DT/tau1;
              St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                    + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
            }
            else
              St = exp( DT/taud)/(p2*p2 * p4*p4) 
                + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
          }
          else 
          {
            if(fabs(p1) < EPS)
            {
              tmp = p3*DT/tau2;
              St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                    + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
            }
            else 
              St = exp(-DT/taud)/(p1*p1 * p3*p3) 
                + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

/*
 * auto-covariance of a line for an array of times with linear amplitude
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
void Sll_array_gamma_linear(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =      pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =      pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau01 - tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
            {
              tmp = p2*DT/tau1;
              St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                    + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
            }
            else
              St = exp( DT/taud)/(p2*p2 * p4*p4) 
                + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
          }
          else 
          {
            if(fabs(p1) < EPS)
            {
              tmp = p3*DT/tau2;
              St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                    + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
            }
            else 
              St = exp(-DT/taud)/(p1*p1 * p3*p3) 
                + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different lines
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_gamma(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau01, tau02, tau1, tau2, taud, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      DT = Dt - (tau01-tau02);
  
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
        {
          tmp = p2*DT/tau1;
          St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
        }
        else
          St = exp( DT/taud)/(p2*p2 * p4*p4) 
            + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
      }
      else 
      {
        if(fabs(p1) < EPS)
        {
          tmp = p3*DT/tau2;
          St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
        }
        else 
          St = exp(-DT/taud)/(p1*p1 * p3*p3) 
            + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

/*
 * covariance between different lines with linear amplitude
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_gamma_linear(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau01, tau02, tau1, tau2, taud, fg12, p1, p2, p3, p4, p51, p52, tmp;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01-tau02);
  
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
        {
          tmp = p2*DT/tau1;
          St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
        }
        else
          St = exp( DT/taud)/(p2*p2 * p4*p4) 
            + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
      }
      else 
      {
        if(fabs(p1) < EPS)
        {
          tmp = p3*DT/tau2;
          St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
        }
        else 
          St = exp(-DT/taud)/(p1*p1 * p3*p3) 
            + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

/*
 * covariance between different lines for an array of times
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
void Sll2_array_gamma(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau01, tau02, tau1, tau2, taud, p1, p2, p3, p4, p51, p52, tmp;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau01-tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
            {
              tmp = p2*DT/tau1;
              St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                    + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
            }
            else
              St = exp( DT/taud)/(p2*p2 * p4*p4) 
                + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
          }
          else 
          {
            if(fabs(p1) < EPS)
            {
              tmp = p3*DT/tau2;
              St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                    + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
            }
            else 
              St = exp(-DT/taud)/(p1*p1 * p3*p3) 
                + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
          }

          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different lines for an array of times with linear amplitude
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
void Sll2_array_gamma_linear(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau01, tau02, tau1, tau2, taud, p1, p2, p3, p4, p51, p52, tmp;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =      pm[idxk1 + 0] ;
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;
    
    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =      pm[idxk2 + 0] ;
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau01-tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
            {
              tmp = p2*DT/tau1;
              St = exp( DT/taud) * ( (tau1/tau2)*(tau1/tau2)/(p2*p2*p2*p2)*(3-2*tmp+tmp*tmp/2) 
                                    + (1-p3*DT/tau2+2*p3/p51)/(p3*p3 * p52*p52) );
            }
            else
              St = exp( DT/taud)/(p2*p2 * p4*p4) 
                + exp( DT/tau2)/(p4*p4 * p52*p52) * (-(1-p4*DT/tau2+2*p4/p51) + (1-p3*DT/tau2+2*p3/p51)*(p4*p4)/(p3*p3));
          }
          else 
          {
            if(fabs(p1) < EPS)
            {
              tmp = p3*DT/tau2;
              St = exp(-DT/taud) * ( (tau2/tau1)*(tau2/tau1)/(p3*p3*p3*p3)*(3+2*tmp+tmp*tmp/2) 
                                    + (1+p2*DT/tau1+2*p2/p52)/(p2*p2 * p51*p51) );
            }
            else 
              St = exp(-DT/taud)/(p1*p1 * p3*p3) 
                + exp(-DT/tau1)/(p1*p1 * p51*p51) * (-(1+p1*DT/tau1+2*p1/p52) + (1+p2*DT/tau1+2*p2/p52)*(p1*p1)/(p2*p2));
          }
          
          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

/*===============================================================================
 * Gamma transfer function with k=1, or exponential transfer function
 * Gamma(x) = A/T * exp[-(x-t0)/T]
 * parameters: A, t0, T
 * constraints: x>t0, T>0
 *===============================================================================
 */
/*
 * covaraince between continuum and line
 *
 */
double Slc_exp(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, Sttot, p1, p2;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);

    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    DT = Dt - tau0;

    if(DT<=0)
    {
      St = exp(DT/taud)/p1;
    }
    else 
    {
      if(fabs(p2) < EPS)
        St = exp(-DT/taud) * (DT/taud + 1.0/p1);
      else
        St = exp(-DT/taud)/p2 + exp(-DT/tau1)*(-1.0/p2 + 1.0/p1);
    }
    
    Sttot += St * fg;
  }

  return Sttot;
}

/*
 * covaraince between continuum and line with linear amplitude
 *
 */
double Slc_exp_linear(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, Sttot, p1, p2;
  int idx, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  Dt = tline - tcon;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);

    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    DT = Dt - tau0;

    if(DT<=0)
    {
      St = exp(DT/taud)/p1;
    }
    else 
    {
      if(fabs(p2) < EPS)
        St = exp(-DT/taud) * (DT/taud + 1.0/p1);
      else
        St = exp(-DT/taud)/p2 + exp(-DT/tau1)*(-1.0/p2 + 1.0/p1);
    }
    
    Sttot += St * fg;
  }

  return Sttot;
}

/*
 * covariance between continuum and line for an array of time
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_exp(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, p1, p2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg = exp(pm[idxk + 0]);
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);
    
    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT<=0)
        {
          St = exp(DT/taud)/p1;
        }
        else 
        {
          if(fabs(p2) < EPS)
            St = exp(-DT/taud) * (DT/taud + 1.0/p1);
          else
            St = exp(-DT/taud)/p2 + exp(-DT/tau1)*(-1.0/p2 + 1.0/p1);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

/*
 * covariance between continuum and line for an array of time with linear amplitude
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array_exp_linear(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double *pm = (double *)model;
  double Dt, DT, taud, fg, tau1, tau0, St, p1, p2;
  int idx, i, j, k, idxk;
  
  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<ncon; i++)
  {
    for(j=0; j<nline; j++)
    {
      Smat[i*nline+j] = 0.0;
    }
  }

  for(k=0; k<num_gaussian; k++)
  {
    idxk = idx + 1 + 3*k;
    fg =     pm[idxk + 0];
    tau0 =   pm[idxk + 1] ;
    tau1 = exp(pm[idxk + 2]);
    
    p1 = (taud + tau1)/taud;
    p2 = (taud - tau1)/taud;

    for(i=0; i<ncon; i++)
    {
      for(j=0; j<nline; j++)
      {
        Dt = tline[j] - tcon[i];

        DT = Dt - tau0;

        if(DT<=0)
        {
          St = exp(DT/taud)/p1;
        }
        else 
        {
          if(fabs(p2) < EPS)
            St = exp(-DT/taud) * (DT/taud + 1.0/p1);
          else
            St = exp(-DT/taud)/p2 + exp(-DT/tau1)*(-1.0/p2 + 1.0/p1);
        }

        Smat[i*nline + j] += St * fg;
      }
    }
  }

  return;
}

/*
 * auto-covariance of a line for an array of time
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
double Sll_exp(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01 - tau02);
    
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
          St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
        else
          St = exp( DT/taud)/(p2 * p4) 
            + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
      }
      else 
      {
        if(fabs(p1) < EPS)
          St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
        else
          St = exp(-DT/taud)/(p1 * p3) 
            + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of a line for an array of time with linear amplitude
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
double Sll_exp_linear(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, Sttot;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52;
  double *pm = (double *)model;
  int idx, k1, k2, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx = idx_line_pm[nds][nls];

  Sttot = 0.0;
  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =      pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =      pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);
      
      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01 - tau02);
    
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
          St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
        else
          St = exp( DT/taud)/(p2 * p4) 
            + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
      }
      else 
      {
        if(fabs(p1) < EPS)
          St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
        else
          St = exp(-DT/taud)/(p1 * p3) 
            + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
      }

      Sttot += St*fg12;
    }
  }
  
  return Sttot;
}

/*
 * auto-covariance of a line for an array of times
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
void Sll_array_exp(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau01 - tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
              St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
            else
              St = exp( DT/taud)/(p2 * p4) 
                + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
          }
          else 
          {
            if(fabs(p1) < EPS)
              St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
            else
              St = exp(-DT/taud)/(p1 * p3) 
                + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

/*
 * auto-covariance of a line for an array of times with linear amplitude
 *
 * nds: index of dataset
 * nls: indexes of line
 *
 */
void Sll_array_exp_linear(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
{
  double Dt, DT, St;
  double taud, fg1, tau01, tau1, fg2, tau02, tau2, fg12, p1, p2, p3, p4, p51, p52;
  double *pm = (double *)model;
  int idx, k1, k2, i, j, idxk1, idxk2;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx = idx_line_pm[nds][nls];

  for(i=0; i<nline; i++)
  {
    for(j=0; j<=i; j++)
    {
      Smat[i*nline + j] = 0.0;
    }
  }

  for(k1=0; k1<num_gaussian; k1++)
  {
    idxk1 = idx + 1 + k1*3;
    fg1 =      pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);
    
    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2=0; k2<num_gaussian; k2++)
    {
      idxk2 = idx + 1 + k2*3;
      fg2 =      pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      for(i=0; i<nline; i++)
      {
        for(j=0; j<=i; j++)
        {
          Dt = tline[i] - tline[j];
          DT = Dt - (tau01 - tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
              St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
            else
              St = exp( DT/taud)/(p2 * p4) 
                + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
          }
          else 
          {
            if(fabs(p1) < EPS)
              St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
            else
              St = exp(-DT/taud)/(p1 * p3) 
                + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
          }

          Smat[i*nline + j] += St * fg12;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different lines
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_exp(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau01, tau02, tau1, tau2, taud, fg12, p1, p2, p3, p4, p51, p52;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;
  
      DT = Dt - (tau01-tau02);
  
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
          St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
        else
          St = exp( DT/taud)/(p2 * p4) 
            + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
      }
      else 
      {
        if(fabs(p1) < EPS)
          St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
        else
          St = exp(-DT/taud)/(p1 * p3) 
            + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

/*
 * covariance between different lines with linear amplitude
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2_exp_linear(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, idxk1, idxk2;
  double fg1, fg2, tau01, tau02, tau1, tau2, taud, fg12, p1, p2, p3, p4, p51, p52;
  double Dt, DT, St, Sttot;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);
  Dt = t1 - t2;

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  Sttot = 0.0;
  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =     pm[idxk1 + 0];
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =     pm[idxk2 + 0];
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      DT = Dt - (tau01-tau02);
  
      if(DT<=0)
      {
        if(fabs(p4) < EPS)
          St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
        else
          St = exp( DT/taud)/(p2 * p4) 
            + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
      }
      else 
      {
        if(fabs(p1) < EPS)
          St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
        else
          St = exp(-DT/taud)/(p1 * p3) 
            + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
      }

      Sttot += St*fg12;
    }
  }

  return Sttot;
}

/*
 * covariance between different lines for an array of times
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
void Sll2_array_exp(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau01, tau02, tau1, tau2, taud, p1, p2, p3, p4, p51, p52;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 = exp(pm[idxk1 + 0]);
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;

    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 = exp(pm[idxk2 + 0]);
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau01-tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
              St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
            else
              St = exp( DT/taud)/(p2 * p4) 
                + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
          }
          else 
          {
            if(fabs(p1) < EPS)
              St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
            else
              St = exp(-DT/taud)/(p1 * p3) 
                + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
          }

          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

/*
 * covariance between different lines for an array of times with linear amplitude
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
void Sll2_array_exp_linear(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
{
  double *pm=(double *)model;
  int idx1, idx2, idx, k1, k2, i, j, idxk1, idxk2;
  double fg1, fg2, fg12, tau01, tau02, tau1, tau2, taud, p1, p2, p3, p4, p51, p52;
  double Dt, DT, St;

  idx = idx_con_pm[nds];
  taud = exp(pm[idx+2]);

  idx1 = idx_line_pm[nds][nls1];
  idx2 = idx_line_pm[nds][nls2];

  for(i=0; i<nline1; i++)
  {
    for(j=0; j<nline2; j++)
    {
      Smat[i*nline2 + j] = 0.0;
    }
  }

  for(k1 = 0; k1<num_gaussian; k1++)
  {
    idxk1 = idx1 + 1 + k1*3;
    fg1 =      pm[idxk1 + 0] ;
    tau01 =    pm[idxk1 + 1] ;
    tau1 = exp(pm[idxk1 + 2]);

    p1 = (taud-tau1)/taud;
    p2 = (taud+tau1)/taud;
    
    for(k2 = 0; k2 < num_gaussian; k2++)
    {
      idxk2 = idx2 + 1 + k2*3;
      fg2 =      pm[idxk2 + 0] ;
      tau02 =    pm[idxk2 + 1] ;
      tau2 = exp(pm[idxk2 + 2]);

      fg12 = fg1 * fg2;
      p3 = (taud+tau2)/taud;
      p4 = (taud-tau2)/taud;
      p51 = (tau1+tau2)/tau1;
      p52 = (tau1+tau2)/tau2;

      for(i=0; i<nline1; i++)
      {
        for(j=0; j<nline2; j++)
        {
          Dt = tline1[i] - tline2[j];
          DT = Dt - (tau01-tau02);

          if(DT<=0)
          {
            if(fabs(p4) < EPS)
              St = exp( DT/taud) * ( -DT/tau2/p2 + tau1/tau2*1.0/p2/p2 + 1.0/(p3*p52) );
            else
              St = exp( DT/taud)/(p2 * p4) 
                + exp( DT/tau2)/(p4 * p52) * (-1.0 + p4/p3);
          }
          else 
          {
            if(fabs(p1) < EPS)
              St = exp(-DT/taud) * ( DT/tau1/p3 + tau2/tau1*1.0/p3/p3 + 1.0/(p2*p51) );
            else
              St = exp(-DT/taud)/(p1 * p3) 
                + exp(-DT/tau1)/(p1 * p51) * (-1.0 + p1/p2);
          }
          
          Smat[i*nline2 + j] += St*fg12;
        }
      }
    }
  }

  return;
}

/*=========================================================*/
void test_covariance()
{
  test_gamma();
  test_gauss();
  test_tophat();
  test_exp();
}

/*
 * plot covariance functions
 */
void test_gamma()
{
  FILE *fp;
  void *model = (void *)malloc(7*sizeof(double));
  double *pm = (double *)model;
  int nt, i;
  double dt, Scc, Slc, Sll, Slc_single, Sll_single, taud, width; 

  num_gaussian = 1;
  idx_con_pm = malloc(1*sizeof(double));
  idx_line_pm = malloc(1*sizeof(double*));
  idx_line_pm[0] = malloc(1*sizeof(double));
  idx_con_pm[0] = 0;
  idx_line_pm[0][0] = 3;

  nt=5000;

  pm[2] = log(20);
  pm[4] = log(1.0);
  pm[5] = 30.0;
  pm[6] = log(20.0);

  taud = exp(pm[2]);
  
  fp = fopen("data/S_gamma.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = -200.0 + 400.0/(nt-1) * i;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_gamma(0.0, dt, model, 0, 0);
    Sll = Sll_gamma(0.0, dt, model, 0, 0);
    Slc_single = Slc_single_gamma(0.0, dt, model, 0, 0, 0);
    Sll_single = Sll_single_gamma(0.0, dt, model, 0, 0, 0);
    fprintf(fp, "%f %e %e %e %e %e\n", dt, Scc, Slc, Sll, Slc_single, Sll_single);
  }
  fclose(fp);

  fp = fopen("data/S_gamma_width.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = 0.0;
    width = log(exp(pm[2]) + (-0.01+0.02/(nt-1)*i));
    pm[6] = width;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_gamma(0.0, dt, model, 0, 0);
    Sll = Sll_gamma(0.0, dt, model, 0, 0);
    Slc_single = Slc_single_gamma(0.0, dt, model, 0, 0, 0);
    Sll_single = Sll_single_gamma(0.0, dt, model, 0, 0, 0);
    fprintf(fp, "%f %e %e %e %e %e\n", width, Scc, Slc, Sll, Slc_single, Sll_single);
  }
  fclose(fp);

  free(model);
  free(idx_con_pm);
  free(idx_line_pm[0]);
  free(idx_line_pm);
}

void test_exp()
{
  FILE *fp;
  void *model = (void *)malloc(7*sizeof(double));
  double *pm = (double *)model;
  int nt, i;
  double dt, Scc, Slc, Sll, Slc_single, Sll_single, taud, width; 

  num_gaussian = 1;
  idx_con_pm = malloc(1*sizeof(double));
  idx_line_pm = malloc(1*sizeof(double*));
  idx_line_pm[0] = malloc(1*sizeof(double));
  idx_con_pm[0] = 0;
  idx_line_pm[0][0] = 3;

  nt=500;

  pm[2] = log(20);
  pm[4] = log(1.0);
  pm[5] = 30.0;
  pm[6] = log(20.0);

  taud = exp(pm[2]);
  
  fp = fopen("data/S_exp.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = -200.0 + 400.0/(nt-1) * i;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_exp(0.0, dt, model, 0, 0);
    Sll = Sll_exp(0.0, dt, model, 0, 0);
    Slc_single = Slc_single_exp(0.0, dt, model, 0, 0, 0);
    Sll_single = Sll_single_exp(0.0, dt, model, 0, 0, 0);
    fprintf(fp, "%f %e %e %e %e %e\n", dt, Scc, Slc, Sll, Slc_single, Sll_single);
  }
  fclose(fp);

  fp = fopen("data/S_exp_width.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = 0.0;
    width = log(exp(pm[2]) + (-0.01+0.02/(nt-1)*i));
    pm[6] = width;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_exp(0.0, dt, model, 0, 0);
    Sll = Sll_exp(0.0, dt, model, 0, 0);
    Slc_single = Slc_single_exp(0.0, dt, model, 0, 0, 0);
    Sll_single = Sll_single_exp(0.0, dt, model, 0, 0, 0);
    fprintf(fp, "%f %e %e %e %e %e\n", width, Scc, Slc, Sll, Slc_single, Sll_single);
  }
  fclose(fp);

  free(model);
  free(idx_con_pm);
  free(idx_line_pm[0]);
  free(idx_line_pm);
  
}

void test_gauss()
{
  FILE *fp;
  void *model = (void *)malloc(7*sizeof(double));
  double *pm = (double *)model;
  int nt, i;
  double dt, Scc, Slc, Sll, taud; 

  num_gaussian = 1;
  idx_con_pm = malloc(1*sizeof(double));
  idx_line_pm = malloc(1*sizeof(double*));
  idx_line_pm[0] = malloc(1*sizeof(double));
  idx_con_pm[0] = 0;
  idx_line_pm[0][0] = 3;

  nt=500;

  pm[2] = log(50);
  pm[4] = log(1.0);
  pm[5] = 30.0;
  pm[6] = log(20.0);

  taud = exp(pm[2]);
  
  fp = fopen("data/S_gauss.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = -200.0 + 400.0/(nt-1) * i;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_gauss(0.0, dt, model, 0, 0);
    Sll = Sll_gauss(0.0, dt, model, 0, 0);
    fprintf(fp, "%f %e %e %e\n", dt, Scc, Slc, Sll);
  }

  free(model);
  free(idx_con_pm);
  free(idx_line_pm[0]);
  free(idx_line_pm);
  fclose(fp);
}

void test_tophat()
{
  FILE *fp;
  void *model = (void *)malloc(7*sizeof(double));
  double *pm = (double *)model;
  int nt, i;
  double dt, Scc, Slc, Sll, taud; 

  num_gaussian = 1;
  idx_con_pm = malloc(1*sizeof(double));
  idx_line_pm = malloc(1*sizeof(double*));
  idx_line_pm[0] = malloc(1*sizeof(double));
  idx_con_pm[0] = 0;
  idx_line_pm[0][0] = 3;

  nt=500;

  pm[2] = log(50);
  pm[4] = log(1.0);
  pm[5] = 30.0;
  pm[6] = log(20.0);

  taud = exp(pm[2]);
  
  fp = fopen("data/S_tophat.txt", "w");
  for(i=0; i<nt; i++)
  {
    dt = -200.0 + 400.0/(nt-1) * i;
    Scc = exp(-fabs(dt)/taud);
    Slc = Slc_tophat(0.0, dt, model, 0, 0);
    Sll = Sll_tophat(0.0, dt, model, 0, 0);
    fprintf(fp, "%f %e %e %e\n", dt, Scc, Slc, Sll);
  }

  free(model);
  free(idx_con_pm);
  free(idx_line_pm[0]);
  free(idx_line_pm);
  fclose(fp);
}
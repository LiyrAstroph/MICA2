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

#include "proto.h"

void *best_model_line;   /*!< best model */
void *best_model_std_line;  /*!< standard deviation of the best model */

void mc_line()
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

  mc_line_init();

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
      printf("# number of Gaussian: %d\n", num_gaussian);
    }
    
    sprintf(postfix, "_%d", num_gaussian);
    strcpy(argv[argc], "-x");
    strcpy(argv[argc+1], postfix);

    logz_arr[j] = dnest_line(argc+2, argv);

    postprocess_line();

    output_reconstrction(best_model_line);

    if(logz_max < logz_arr[j])
    {
      logz_max = logz_arr[j];
      jzmax = j;
    }
  }
  
  if(thistask == roottask)
  {
    printf("*****************************************************\n");
    for(j=0; j<parset.num_gaussian_diff; j++)
    {
      printf("number of Gaussian: %d, evidence: %f\n", parset.num_gaussian_low + j, logz_arr[j]);
    }

    printf("best number of Gaussian: %d.\n", parset.num_gaussian_low + jzmax);
    printf("*****************************************************\n");
  }  

  mc_line_end();

  /* clear up argv */
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

/* 
 * reconstruct line lcs acording to the time span of the continuum.
 */
void output_reconstrction(const void *model)
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

      nall[i][0] = (int) fmin(dataset[i].con.n*nscale, 500);
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

    for(i=0; i<nset; i++)
    {
      tbeg = dataset[i].con.t[0];
      tend = dataset[i].con.t[dataset[i].con.n-1];
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
        recostruct_line_from_varmodel(post_model, i, nall[i], tall[i], fall[i], feall[i], yq); 

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

      fprintf(fp, "# %d",nall[i][0]);
      for(j=0;  j < dataset[i].nlset; j++)
        fprintf(fp, ":%d", nall[i][1+j]);
      fprintf(fp, "\n");
      
      /* output reconstructed continuum */
      for(k=0; k<nall[i][0]; k++)
        fprintf(fp, "%f %f %f\n", tall[i][k], fall_best[i][k] * dataset[i].con.scale, fall_std[i][k] * dataset[i].con.scale);
      fprintf(fp, "\n");

      /* output reconstructed lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
          fprintf(fp, "%f %f %f\n", tall[i][np+k], fall_best[i][np+k] * dataset[i].line[j].scale, fall_std[i][np+k] * dataset[i].line[j].scale);
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

/* 
 * reconstruct line lcs acording to their own time spans.
 */
void output_reconstrction2(const void *model)
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
    void *post_model;
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

    for(i=0; i<nset; i++)
    {
      /* time nodes of continuum */
      nall[i][0] = (int)fmin(dataset[i].con.n * nscale, 500);
      tspan = dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0];
      for(j=0; j<nall[i][0]; j++)
        tall[i][j] = (tspan+0.1*tspan)/(nall[i][0]-1.0) * j 
                        + dataset[i].con.t[0]-0.05*tspan;

      /* time nodes of lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        nall[i][1+j] = (int)fmin(dataset[i].line[j].n * nscale, 500);
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
      
      for(i=0; i<nset; i++)
      {
        /* reconstuct all the light curves */
        recostruct_line_from_varmodel(post_model, i, nall[i], tall[i], fall[i], feall[i], yq); 

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
        fprintf(fp, "%f %f %f\n", tall[i][k], fall_best[i][k] * dataset[i].con.scale, fall_std[i][k] * dataset[i].con.scale);
      fprintf(fp, "\n");

      /* output reconstructed lines */
      np = nall[i][0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[i][1+j]; k++)
          fprintf(fp, "%f %f %f\n", tall[i][np+k], fall_best[i][np+k] * dataset[i].line[j].scale, fall_std[i][np+k] * dataset[i].line[j].scale);
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
  
  best_model_line = malloc(size_of_modeltype);
  best_model_std_line = malloc(size_of_modeltype);
  
  if(thistask == roottask)
  {
    FILE *fp;

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
  }
  return;
}

/*!
 * matrix operation A^-1 x B is implemented by calling functions
 *    inverse_mat()
 *    multiply_mat_MN()
 */
void recostruct_line_from_varmodel2(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall)
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
  y = ybuf + nall_data;
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
void recostruct_line_from_varmodel(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq, *yave;
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
  y = ybuf + nall_data;
  Cq = y + nall_data;
  yq = Cq + nqall*nqall;
  yave = yq + nqall;
  yuq = yave + nall_data; 
  Larr_rec = yuq + ntall;

  PEmat1 = malloc(ntall * ntall * sizeof(double));
  PEmat2 = malloc(ntall * ntall * sizeof(double));
  PEmat3 = malloc(ntall * ntall * sizeof(double));
  PEmat4 = malloc(ntall * ntall * sizeof(double));

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
      if(flag_warning)printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
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
      if(flag_warning)printf("prob >0!\n");
      return prob;
    }
    
    prob += prob1 -0.5*lndet - 0.5*lndet_ICq;
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

  free(best_model_line);
  free(best_model_std_line);
  free(Smat_lc);
  free(Smat_ll);
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

/*
 * auto-covariance of line
 *
 * nds: index of dataset
 * nls: index of line
 *
 */
double Sll(double t1, double t2, const void *model, int nds, int nls)
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


/*
 * auto-covariance of line for an array of times
 *
 * nds: index of dataset
 * nls: index of line
 *
 */

void Sll_array(double *tline, int nline, const void *model, int nds, int nls, double *Smat)
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

/*
 * covariance between different line
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */
double Sll2(double t1, double t2, const void *model, int nds, int nls1, int nls2)
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

/*
 * covariance between different line for an array of times
 *
 * nds: index of dataset
 * nls1, nls2: indexes of line1 and line2
 *
 */

void Sll2_array(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat)
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

/*
 * covariance between continuum and line 
 *
 * nds: dataset index
 * nls: line set index
 */
double Slc(double tcon, double tline, const void *model, int nds, int nls)
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

/*
 * covariance between continuum and line for an array of times
 *
 * nds: dataset index
 * nls: line set index
 */
void Slc_array(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat)
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


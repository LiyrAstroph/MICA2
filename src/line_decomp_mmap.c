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

#include "dnest.h"
#include "dnest_line.h"
#include "allvars.h"

#include "proto.h"

void output_decompose_line_parallel_mmap()
{
  FILE *fp=NULL, *fp_sample=NULL;
  char fname[200];
  int i, j, k, m, kgau;
  double **tall, **fall, **feall, **feall_max, **fall_best, **fall_std, *yq;
  double **fall_best_buf, **fall_std_buf, **feall_max_buf;
  int **nall, *ntall, np;
    
  int num_ps, size_of_modeltype, num_ps_task, num_ps_task_average;
  void *posterior_sample_all, *posterior_sample_task;
  double *post_model, *ps;
  char posterior_sample_file[MICA_MAX_STR_LENGTH];
  
  /* time span of reconstruction */
  double tspan;
  double tbeg, tend;

  size_of_modeltype = num_params * sizeof(double);
  
  if(thistask == roottask)
  {
    printf("start mmap decompose...\n");
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

    /* note here read all posterior sample */
    posterior_sample_all = (void *)malloc(num_ps * size_of_modeltype);
    post_model = (double *)posterior_sample_all;

    for(m=0; m<num_ps; m++)
    {      
      // read sample
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp_sample, "%lf", post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp_sample, "\n");
      
      /* in pmap, need to transform response ratio into normal values as in gmodel */
      if(parset.model == pmap)
      {
        transform_response_ratio_inplace(post_model);
      }
      
      post_model += num_params;
    }

    fclose(fp_sample);
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
 
  nall = malloc(nset*sizeof(int *));
  tall = malloc(nset*sizeof(double *));
  fall = malloc(nset*sizeof(double *));
  feall = malloc(nset*sizeof(double *));
  feall_max = malloc(nset*sizeof(double *));
  fall_best = malloc(nset*sizeof(double *));
  fall_std = malloc(nset*sizeof(double *));
  ntall = malloc(nset*sizeof(int));
    
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
    tall[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall[i] = (double *)malloc(ntall[i] * sizeof(double));
    feall[i] = (double *)malloc(ntall[i] * sizeof(double));
    feall_max[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_best[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_std[i] = (double *)malloc(ntall[i] * sizeof(double));
  }
  
  yq = (double *)malloc(nq*(1+nlset_max)*sizeof(double));
  
  /* bufs for reduce */
  feall_max_buf = (double **)malloc(nset*sizeof(double *));
  fall_best_buf = (double **)malloc(nset*sizeof(double *));
  fall_std_buf = (double **)malloc(nset*sizeof(double *));
  for(i=0; i<nset; i++)
  {
    feall_max_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_best_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
    fall_std_buf[i] = (double *)malloc(ntall[i] * sizeof(double));
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
 
  for(kgau=0; kgau<num_gaussian; kgau++)
  {
    if(thistask == roottask )
    {
      printf("%d-th component\n", kgau);
      pb_init(pb_mica, '#', 50, num_ps_task);
      showPercent(pb_mica, true);
      showCount(pb_mica, false);
      pb_update(pb_mica, 0);
      pb_print(pb_mica);

      sprintf(fname, "%s/%s%s_comp%d", parset.file_dir, "data/pline.txt", postfix, kgau);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      fprintf(fp, "# %d\n", nset);
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
    }
      
    ps = (double *)posterior_sample_task;
    for(m=0; m<num_ps_task; m++)
    {
      if(thistask == roottask)
      {
        pb_update(pb_mica, m+1);
        pb_print(pb_mica);
      }
      
      // printf("# %d-th component of sample %d on task %d\n", kgau+1, m, thistask);      
      for(i=0; i<nset; i++)
      {
        /* reconstuct all the light curves */
        decompose_single_component_mmap((void *)ps, i, nall[i], tall[i], fall[i], feall[i], yq, kgau); 

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
      }
      ps += num_params;
    }

    /* reduce the reconstruction */
    for(i=0; i<nset; i++)
    {
      MPI_Reduce(fall_best[i], fall_best_buf[i], ntall[i], MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
      MPI_Reduce(fall_std[i], fall_std_buf[i], ntall[i], MPI_DOUBLE, MPI_SUM, roottask, MPI_COMM_WORLD);
      MPI_Reduce(feall_max[i], feall_max_buf[i], ntall[i], MPI_DOUBLE, MPI_MAX, roottask, MPI_COMM_WORLD);
    }

    if(thistask == roottask)
    {
      for(i=0; i<nset; i++)
      { 
        for(k=0; k<nall[i][0]; k++)
        {
          fall_best_buf[i][k] /= num_ps;
          fall_std_buf[i][k] /= num_ps;
          fall_std_buf[i][k] = sqrt(fall_std_buf[i][k] - fall_best_buf[i][k]*fall_best_buf[i][k]);
          fall_std_buf[i][k] = fmax(fall_std_buf[i][k], feall_max_buf[i][k]);
        }
          
        np = nall[i][0];
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nall[i][1+j]; k++)
          {
            fall_best_buf[i][np+k] /= num_ps;
            fall_std_buf[i][np+k] /= num_ps;
            fall_std_buf[i][np+k] = sqrt(fall_std_buf[i][np+k] - fall_best_buf[i][np+k]*fall_best_buf[i][np+k]);
            fall_std_buf[i][np+k] = fmax(fall_std_buf[i][np+k], feall_max_buf[i][np+k]);
          }          
          np += nall[i][1+j];
        }

        fprintf(fp, "# %d",nall[i][0]);
        for(j=0;  j < dataset[i].nlset; j++)
          fprintf(fp, ":%d", nall[i][1+j]);
        fprintf(fp, "\n");
        
        /* output reconstructed continuum */
        for(k=0; k<nall[i][0]; k++)
          fprintf(fp, "%f %f %f\n", tall[i][k], fall_best_buf[i][k] * dataset[i].con.scale, fall_std_buf[i][k] * dataset[i].con.scale);
        fprintf(fp, "\n");

        /* output reconstructed lines */
        np = nall[i][0];
        for(j=0; j<dataset[i].nlset; j++)
        {
          for(k=0; k<nall[i][1+j]; k++)
            fprintf(fp, "%f %f %f\n", tall[i][np+k], fall_best_buf[i][np+k] * dataset[i].line[j].scale, fall_std_buf[i][np+k] * dataset[i].line[j].scale);
          fprintf(fp, "\n");
          np += nall[i][1+j];
        }  
      }
      fclose(fp);
    }
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
    free(fall_best[i]);
    free(fall_std[i]);

    free(feall_max_buf[i]);
    free(fall_best_buf[i]);
    free(fall_std_buf[i]);
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

  free(feall_max_buf);
  free(fall_best_buf);
  free(fall_std_buf);

  if(thistask == roottask)
  {
    printf("end mmap decompose...\n");
  }

  return;
}

/* 
 * decompose line component from each Gaussian
 */
void decompose_single_component_mmap(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall, int kgau)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq, *yave;
  int i, j, k, m, info, idx, *ipiv;
  double *PEmat1, *PEmat2, *PEmat3, *PEmat4;
  int nall_data, nqall, ntall, np, ntall_max;
  double *fall_data;
  double sigma, tau, *pm=(double *)model;
  double fq_weight, *fq_all; /* q weight of each component */

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

  /* calculate total responses of each dataset */
  fq_all = malloc((dataset[nds].nlset)*sizeof(double));
  for(j=0; j<dataset[nds].nlset; j++)
  {
    idx = idx_line_pm[nds][j];
    fq_all[j] = 0.0;
    if(parset.flag_negative_resp == 0)
    {
      for(k=0; k<num_gaussian; k++)
      {
        fq_all[j] += exp(pm[idx+1+k*3]);
      }
    }
  }
  
  /* determine the maximum size */
  ntall_max = ntall>nall_max?ntall:nall_max;
  
  PEmat1 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat2 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat3 = malloc(ntall_max * ntall_max * sizeof(double));
  PEmat4 = malloc(ntall_max * ntall_max * sizeof(double));

  set_covar_Pmat_data_line_array_mmap(model, nds);
  
  set_covar_Umat_line_single_mmap(model, nds, nall, tall, kgau);
  
  set_covar_Amat_line_single_mmap(model, nds, nall, tall, kgau);

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
    idx = idx_line_pm[nds][k];
    if(parset.flag_negative_resp == 0)
    {
      fq_weight = exp(pm[idx+1+3*kgau])/fq_all[k];
      for(i=0; i<nall[1+k]; i++)
      {
        for(j=0; j<nqall; j++)
        {
          Larr_rec[(np+i)*nqall + j] = 0.0;
        }

        Larr_rec[(np+i)*nqall + nq + k*nq + 0] = 1.0 * fq_weight;
        for(j=1; j<nq; j++)
          Larr_rec[(np+i)*nqall + nq + k*nq + j] = pow(tall[np+i], j) * fq_weight;
      }
    }
    else 
    {
      /* do not add long-term trend when there is negative response */
      for(i=0; i<nall[1+k]; i++)
      {
        for(j=0; j<nqall; j++)
        {
          Larr_rec[(np+i)*nqall + j] = 0.0;
        }
      }
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
  //multiply_mat_MN_transposeA(PEmat2, Larr, PEmat3, ntall, nqall, nall_data);
  /* S x C^-1 x L - L */
  //for(i=0; i<ntall*nqall; i++)PEmat3[i] -= Larr_rec[i];
  
  //inverse_mat(Cq, nqall, &info, ipiv);

  //multiply_mat_MN(PEmat3, Cq, PEmat2, ntall, nqall, nqall);
  /* (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  //multiply_mat_MN_transposeB(PEmat2, PEmat3, PEmat4, ntall, ntall, nqall);
  
  for(i=0; i<ntall; i++)
  {
    //feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat1[i*ntall+i] + PEmat4[i*ntall + i]);
    feall[i] = sigma * sqrt(ASmat[i*ntall + i] - PEmat1[i*ntall+i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  free(fq_all);
  return;
}

/*!
 * this function sets the covariance matrix at time of data points and reconstruction points
 * for a given component
 */
void set_covar_Umat_line_single_mmap(const void *model, int nds, int *nall, double *tall, int kgau)
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

    /* lines, for continuum (data) - line (rec) covariance, all line components should be included */
    np = dataset[nds].con.n;
    for(k=0; k<dataset[nds].nlset; k++)
    {
      for(j=0; j < dataset[nds].line[k].n; j++)
      {
        t2 = dataset[nds].line[k].t[j];
        USmat[ i * nall_data + j + np ] = Slc_mmap(t1, t2, model, nds, k);
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

      /* continum, for line (data) - continuum (rec), only the component included */
      for(j=0; j<dataset[nds].con.n; j++)
      {
        t2 = dataset[nds].con.t[j];
        USmat[ (i+npline) * nall_data + j ] = Slc_single_mmap(t2, t1, model, nds, k, kgau);
      }

      /* line different */
      np = dataset[nds].con.n;
      for(m=0; m<k; m++)
      {
        for(j=0; j < dataset[nds].line[m].n; j++)
        {
          t2 = dataset[nds].line[m].t[j];
          USmat[ (i + npline) * nall_data + j + np ] = Sll2_single_mmap(t1, t2, model, nds, k, m, kgau);
        }
        np += dataset[nds].line[m].n;
      }
      
      /* line self */
      for(j=0; j<dataset[nds].line[k].n; j++)
      {
        t2 = dataset[nds].line[k].t[j];
        USmat[ (i+npline) * nall_data + j + np ] = Sll_single_mmap(t1, t2, model, nds, k, kgau);
      }
      np += dataset[nds].line[k].n;

      /* line different */
      for(m=k+1; m<dataset[nds].nlset; m++)
      {
        for(j=0; j < dataset[nds].line[m].n; j++)
        {
          t2 = dataset[nds].line[m].t[j];
          USmat[ (i + npline) * nall_data + j + np ] = Sll2_single_mmap(t1, t2, model, nds, k, m, kgau);
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
 * for a given component
 */
void set_covar_Amat_line_single_mmap(const void *model, int nds, int *nall, double *tall, int kgau)
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
        ASmat[i * ntall + np+j] = ASmat[(np+j)*ntall + i] = Slc_mmap(t1, t2, model, nds, k);
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
        ASmat[(j+np)*ntall + np+i] = ASmat[(i+np)*ntall + np+j] = Sll_single_mmap(t1, t2, model, nds, k, kgau);
      }
      ASmat[(j+np)*ntall + (j+np)] = Sll_single2(t1, t1, model, nds, k, kgau);
      
      /* line different */
      npline = np + nall[1+k];
      for(m = k+1; m<dataset[nds].nlset; m++)
      {
        for(i=0; i<nall[1+m]; i++)
        {
          t2 = tall[npline + i];
          ASmat[(j+np)*ntall + npline + i] = ASmat[ (npline + i)*ntall + j+np ] = Sll2_single_mmap(t1, t2, model, nds, k, m, kgau) ;
        }
        npline += nall[1+m];
      }
    }
    np += nall[1+k];
  }
  return;
}

double Slc_single_mmap(double tcon, double tline, const void *model, int nds, int nls, int k_comp)
{
  int i, k, type;
  double Stot;
  
  Stot = 0.0;
  k = k_comp;
  {
    type = parset.str_type_tf_mix[k] - '0';

    Stot += FP_Slc[type](tcon, tline, model, nds, nls, k);
  }
  return Stot;
}

double Sll_single_mmap(double t1, double t2, const void *model, int nds, int nls, int k_comp)
{
  int k1, k2, type1, type2;
  double Stot;

  Stot = 0.0;
  k1 = k_comp;
  {
    type1 = parset.str_type_tf_mix[k1] - '0'; /* type of k1 component */ 
    
    for(k2=0; k2<num_gaussian; k2++)
    {
      type2 = parset.str_type_tf_mix[k2] - '0'; /* type of k2 component */ 

      Stot += FP_Sll[type1*4 + type2](t1, t2, model, nds, nls, nls, k1, k2);
    }
  }

  return Stot;
}

double Sll2_single_mmap(double t1, double t2, const void *model, int nds, int nls1, int nls2, int k_comp)
{
  int k1, k2, type1, type2;
  double Stot;

  Stot = 0.0;
  k1 = k_comp;
  {
    type1 = parset.str_type_tf_mix[k1] - '0'; /* type of k1 component */

    for(k2=0; k2<num_gaussian; k2++)
    {
      type2 = parset.str_type_tf_mix[k2] - '0'; /* type of k2 component */ 

      Stot += FP_Sll[type1*4 + type2](t1, t2, model, nds, nls1, nls2, k1, k2);
    }
  }
  return Stot;
}
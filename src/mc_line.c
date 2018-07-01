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
#include <gsl/gsl_integration.h>

#include "dnest_line.h"
#include "allvars.h"

#include "proto.h"

void *best_model_line;   /*!< best model */
void *best_model_std_line;  /*!< standard deviation of the best model */

void mc_line()
{
  int i, argc=0;
  char **argv;

  argv = malloc(9*sizeof(char *));
  for(i=0; i<9; i++)
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

  strcpy(argv[argc++], "-p");

  mc_line_init();

  dnest_line(argc, argv);

  postprocess_line();

  if(thistask == roottask)
  {
    FILE *fp;
    char fname[200];
    int i, j, k;
    double *tall, *fall, *feall;
    int *nall, ntall, np;

    sprintf(fname, "%s/%s", parset.file_dir, "data/pall.txt");
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
 
    nall = malloc((1+nlset_max) * sizeof(double));
    tall = malloc(nall_max*5 * sizeof(double));
    fall = malloc(nall_max*5 * sizeof(double));
    feall = malloc(nall_max*5 * sizeof(double));

    for(i=0; i<nset; i++)
    {

      /* time nodes of continuum */
      nall[0] = dataset[i].con.n * 5;
      for(j=0; j<nall[0]; j++)
        tall[j] = (dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0]+40.0)/(nall[0]-1.0) * j 
                        + dataset[i].con.t[0]-20.0;

      /* time nodes of lines */
      np = nall[0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        nall[1+j] = dataset[i].line[j].n * 5;
        for(k=0; k<nall[1+j]; k++)
        {
          tall[np+k] = (dataset[i].line[j].t[dataset[i].line[j].n-1] - dataset[i].line[j].t[0]+40.0)/(nall[1+j]-1.0) * k 
                         + dataset[i].line[j].t[0]-20.0;
        }
        np += nall[1+j];
      }

      /* compute total number of points */
      ntall = nall[0];
      for(j=0; j<dataset[i].nlset; j++)
        ntall += nall[j];

      /* reconstuct all the light curves */
      recostruct_line_from_varmodel(best_model_line, i, nall, tall, fall, feall);

      /* output reconstructed continuum */
      for(k=0; k<nall[0]; k++)
        fprintf(fp, "%f %f %f\n", tall[k], fall[k] * dataset[i].con.scale, feall[k] * dataset[i].con.scale);
      fprintf(fp, "\n");

      /* output reconstructed lines */
      np = nall[0];
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<nall[1+j]; k++)
          fprintf(fp, "%f %f %f\n", tall[np+k], fall[np+k] * dataset[i].line[j].scale, feall[np+k] * dataset[i].line[j].scale);
        fprintf(fp, "\n");
        np += nall[1+j];
      }
    }
    
    fclose(fp);

    free(tall);
    free(fall);
    free(feall);
    free(nall);
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


/*!
 *  this function does postprocess for continuum. 
 */
void postprocess_line()
{
  char posterior_sample_file[MICA_MAX_STR_LENGTH];
  double *pm, *pmstd;
  int num_ps, i, j;
  void *posterior_sample, *post_model;
  
  best_model_line = malloc(size_of_modeltype);
  best_model_std_line = malloc(size_of_modeltype);
  
  if(thistask == roottask)
  {
    char fname[200];
    FILE *fp;

    /* get file name of posterior sample file */
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

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

void recostruct_line_from_varmodel(const void *model, int nds, int *nall, double *tall, double *fall, double *feall)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq;
  int i, j, k, m, info;
  double *PEmat1, *PEmat2, *PEmat3, *PEmat4;
  int nall_data, nqall, ntall, np;
  double *tall_data, *fall_data, *feall_data;

  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  nqall = nq * (1  + dataset[nds].nlset);
  nall_data = alldata[nds].n;
  tall_data = alldata[nds].t;
  fall_data = alldata[nds].f;
  feall_data = alldata[nds].fe;

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

  inverse_mat(PCmat, nall_data, &info);

  for(i=0;i<dataset[nds].con.n;i++)
  {
    Larr[i*nqall] = 1.0; 
    for(j=0; j<dataset[nds].nlset; j++)
      Larr[i*nqall + 1 + j] = 0.0;
  }
  np = dataset[nds].con.n;
  for(j=0; j<dataset[nds].nlset; j++)
  {
    for(m=0; m<dataset[nds].line[j].n; m++)
    {
      for(i=0; i<nqall; i++)
        Larr[(np+m)*nqall + i ]  = 0.0;
        
      Larr[(np + m) * nqall + 1 + j] = 1.0;
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
  inverse_mat(Cq, nqall, &info);
  multiply_mat_MN(Cq, yuq, yq, nqall, 1, nqall);

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
    for(j=1; j<nqall; j++)
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
      Larr_rec[(np+i)*nqall + 1+k ] = 1.0;
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
    feall[i] = sqrt(ASmat[i*ntall + i] - PEmat2[i*ntall+i] + PEmat4[i*ntall + i]);
    //printf("%d %f %f %f\n", i, ASmat[i*ntall + i],  PEmat2[i*ntall+i], PEmat4[i*ntall + i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  return;
}

double prob_line_variability(const void *model)
{
  double prob = 0.0;
  int i, j, k, m, np, info, sign;
  double lndet, lndet_ICq;
  double *Larr, *ybuf, *y, *yq, *Cq, *ICq;
  double *tall, *fall, *feall;
  int nall, nqall;

  Larr = workspace;
  ybuf = Larr + nall_max * ((1+nlset_max)*nq);
  y = ybuf + nall_max;
  yq = y + nall_max;
  Cq = yq + (1+nlset_max)*nq;
  ICq = Cq + (1+nlset_max)*nq * (1+nlset_max)*nq;

  /* iterate over all datasets */
  for(k=0; k<2; k++)
  {
    nall = alldata[k].n;
    tall = alldata[k].t;
    fall = alldata[k].f;
    feall = alldata[k].fe;
    nqall = nq * (1+dataset[k].nlset);

    for(i=0;i<dataset[k].con.n;i++)
    {
      Larr[i*nqall] = 1.0; 
      for(j=0; j<dataset[k].nlset; j++)
        Larr[i*nqall + 1 + j] = 0.0;
    }
    np = dataset[k].con.n;
    for(j=0; j<dataset[k].nlset; j++)
    {
      for(m=0; m<dataset[k].line[j].n; m++)
      {
        for(i=0; i<nqall; i++)
          Larr[(np+m)*nqall + i ]  = 0.0;
        
        Larr[(np + m) * nqall + 1 + j] = 1.0;
      }
      np += dataset[k].line[j].n;
    }
     
    set_covar_Pmat_data_line(model, k);
    memcpy(IPCmat, PCmat, nall*nall*sizeof(double));

    inverse_mat(IPCmat, nall, &info); /* calculate C^-1 */
    
    /* calculate L^T*C^-1*L */
    multiply_mat_MN(IPCmat, Larr, ybuf, nall, nqall, nall);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nqall, nqall, nall);
    memcpy(ICq, Cq, nqall*nqall*sizeof(double));

    /* calculate L^T*C^-1*y */
    multiply_matvec(IPCmat, fall, nall, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nqall, 1, nall);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_mat(Cq, nqall, &info);
    multiply_mat_MN(Cq, yq, ybuf, nqall, 1, nqall);
  
    multiply_matvec_MN(Larr, nall, nqall, ybuf, y);
    for(i=0; i<nall; i++)
    {
      y[i] = fall[i] - y[i];
    }
  
    /* y^T x C^-1 x y */
    multiply_matvec(IPCmat, y, nall, ybuf);
    prob += -0.5 * cblas_ddot(nall, y, 1, ybuf, 1);
    
    lndet = lndet_mat3(PCmat, nall, &info, &sign);
    if(info!=0|| sign==-1)
    {
      prob = -DBL_MAX;
      printf("lndet_C %f %d!\n", lndet, sign);
      return prob;
    }
    lndet_ICq = lndet_mat3(ICq, nqall, &info, &sign);
    if(info!=0 || sign==-1 )
    {
      prob = -DBL_MAX;
      printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
      return prob;
    }
    prob += -0.5*lndet - 0.5*lndet_ICq;
  }

  return prob;
}

int mc_line_init()
{
  int i, j, np;

  sprintf(dnest_options_file, "%s/src/OPTIONS1D", parset.file_dir);

  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  perturb_accept = malloc(parset.num_particles * sizeof(int));
  for(i=0; i<parset.num_particles; i++)
  {
    perturb_accept[i] = 0;
  }

  idx_line_pm = malloc(nset * sizeof(double *));
  for(i=0; i<nset; i++)
  {
    idx_line_pm[i] = malloc(nlset_max * sizeof(int));
  }

  if(parset.flag_uniform_tranfuns == 1)
  {
    for(i=0; i<nset; i++)
    {
      for(j=0; j<dataset[i].nlset; j++)
      {
        idx_line_pm[i][j] = 4*j + num_params_var ;
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
        idx_line_pm[i][j] = 4*np + num_params_var;
        //printf("%d %d %d\n", i, j, idx_line_pm[i][j]);
        np++;
      }
    }
    //exit(0);
  }
  return 0;
}

int mc_line_end()
{
  int i;

  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);
  
  free(perturb_accept);

  for(i=0; i<nset; i++)
  {
    free(idx_line_pm[i]);
  }
  free(idx_line_pm);

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
  int i, j, m, l, nall, np, npline;
  double *pm = (double *)model;
  
  /* total number of point of this dataset */
  nall = alldata[k].n;

  /* set variability parameters */
  if(parset.flag_uniform_var_params == 1)
  {
    syserr = (exp(pm[0])-1.0)*dataset[k].con.error_mean;
    tau = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(tau);
  }
  else
  {
    syserr = (exp(pm[3*k])-1.0)*dataset[k].con.error_mean;
    tau = exp(pm[3*k+2]);
    sigma = exp(pm[3*k+1]) * sqrt(tau);
  }

  /* continuum - continuum  */
  for(i=0; i<dataset[k].con.n; i++)
  { 
    t1 = dataset[k].con.t[i];

    for(j=0; j<i; j++)
    {
      t2 = dataset[k].con.t[j];
      PCmat[i*nall+j] = PCmat[j*nall+i] = 
      PSmat[i*nall+j] = PSmat[j*nall+i] = sigma*sigma* exp (- fabs(t1-t2) / tau );
    }

    PSmat[i*nall+i] = sigma * sigma;
    error = dataset[k].con.fe[i]*dataset[k].con.fe[i] + syserr*syserr;
    PCmat[i*nall+i] = PSmat[i*nall+i] + error;
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

    if(parset.flag_uniform_tranfuns == 1)
    {
      syserr_line = (exp(pm[num_params_var + j * 4]) - 1.0) * dataset[k].line[j].error_mean;
    }
    else
    {
      syserr_line = (exp(pm[num_params_var + (k*nlset_max + j) *4 ]) - 1.0) * dataset[k].line[j].error_mean;
    }

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
      PCmat[(np + i)*nall+(np +i)] = PSmat[(np + i)*nall+(np +i)] + error;

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
void set_covar_Umat_line(const void *model, int nds, int *nall, double *tall)
{
  double sigma, taud, t1, t2;
  int i, j, k, m, ntall, nall_data, np, npline;
  double *pm = (double *)model;

  nall_data = alldata[nds].n;
  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  if(parset.flag_uniform_var_params == 1)
  {
    taud = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(taud);
  }
  else
  {
    taud = exp(pm[3*nds+2]);
    sigma = exp(pm[3*nds+1]) * sqrt(taud);
  }

  /* continuum - continuum/line */
  for(i=0; i<nall[0]; i++)
  {
    t1 = tall[i];
   
    /* continuum self */
    for(j=0; j<dataset[nds].con.n; j++)
    {
      t2 = dataset[nds].con.t[j];
      USmat[i * nall_data + j] = sigma*sigma * exp( - fabs(t1-t2)/taud);
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
  double t1, t2, sigma, taud;
  int i, j, k, m, ntall, np, npline;
  double *pm=(double *)model;
  
  ntall = nall[0];
  for(i=0; i<dataset[nds].nlset; i++)
    ntall += nall[i+1];

  if(parset.flag_uniform_var_params == 1)
  {
    taud = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(taud);
  }
  else
  {
    taud = exp(pm[3*nds+2]);
    sigma = exp(pm[3*nds+1]) * sqrt(taud);
  }
  

  /* continuum - continuum/lines*/
  for(i=0; i<nall[0]; i++)
  {
    t1 = tall[i];

    /* continuum */
    for(j=0; j<i; j++)
    {
      t2 = tall[j];
      ASmat[i * ntall + j] = ASmat[j * ntall + i] = sigma*sigma * exp( - fabs(t1-t2)/taud);
    }
    ASmat[i * ntall + i] = sigma * sigma;
    
    /* lines */
    np = nall[0];
    for(k=0; k<dataset[nds].nlset; k++)
    {
      for(j=0; j<nall[1+k]; j++)
      {
        t2 = tall[np + j];
        ASmat[i * ntall + np+j] = ASmat[(np+j)*ntall + i] = Slc(t1, t2, model, nds, k);
      }
      np += nall[k];
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


double Sll(double t1, double t2, const void *model, int nds, int nls)
{
  double Dt, DT, St, lnSt, x1=0.0, x2=0.0;
  double sigma, taud, fg, tau0, wg;
  double *pm = (double *)model;
  int idx;


  if(parset.flag_uniform_var_params == 1)
  {
    taud = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(taud);
  }
  else
  {
    taud = exp(pm[3*nds+2]);
    sigma = exp(pm[3*nds+1]) * sqrt(taud);
  }

  idx = idx_line_pm[nds][nls];
  fg = exp(pm[idx + 1]);
  tau0 =   pm[idx + 2] ;
  wg = exp(pm[idx + 3]);

  Dt = t1 - t2;
  DT = Dt;
  /* 
   * for very large x 
   * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
   * 
   * (see https://en.wikipedia.org/wiki/Error_function)
   *
   * with this equation, the factor exp(wg*wg/taud/taud) will be canceled out.
   */
  if( wg/taud >= 20.0 )
  {
    x1 = -DT/2.0/wg + wg/taud;
    x2 =  DT/2.0/wg + wg/taud;

    St = exp( - DT*DT/wg/wg/4.0 ) / sqrt(PI) * (1.0/x1 * ( 1.0 - 1.0/(2*x1*x1) +  3.0/pow(2*x1*x1, 2) - 15.0/pow(2*x1*x1, 3)) 
                                               +1.0/x2 * ( 1.0 - 1.0/(2*x2*x2) +  3.0/pow(2*x2*x2, 2) - 15.0/pow(2*x2*x2, 3)));
  }
  else
  {
    St = exp(wg*wg/taud/taud) * ( exp(-DT/taud) * erfc( -DT/2.0/wg + wg/taud )
                                 +exp( DT/taud) * erfc(  DT/2.0/wg + wg/taud ) );
  }
  
  St *= 1.0/2.0 * fg*fg * sigma*sigma;
  
  /*if(isnan(St))
  {
    printf("Sll: %f %f %f %f %f %e %f\n", St, sigma, fg, x1, x2, exp(30*30), wg/taud);
    exit(0);
  }*/
  return St;
}

double Sll2(double t1, double t2, const void *model, int nds, int nls1, int nls2)
{
  printf("FFFFF\n");
  exit(0);
  return 0.0;
}
/*
 * nds: dataset index
 * nls: line set index
 */
double Slc(double tcon, double tline, const void *model, int nds, int nls)
{
  double *pm = (double *)model;
  double Dt, DT, sigma, taud, fg, wg, tau0, St, lnSt, x1, x2;
  int idx;
  
  if(parset.flag_uniform_var_params == 1)
  {
    taud = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(taud);
  }
  else
  {
    taud = exp(pm[3*nds+2]);
    sigma = exp(pm[3*nds+1]) * sqrt(taud);
  }

  idx = idx_line_pm[nds][nls];
  fg = exp(pm[idx + 1]);
  tau0 =   pm[idx + 2] ;
  wg = exp(pm[idx + 3]);

  Dt = tline - tcon;
  DT = Dt - tau0;

  /* 
   * for very large x 
   * erfc(x) = exp(-x^2)/(x * pi^1/2) * ( 1 - 1/(2*x*x) + 3/(2*x*x)^2 ... )
   * 
   * (see https://en.wikipedia.org/wiki/Error_function)
   *
   * with this equation, the factor exp(wg*wg/taud/taud/2) will be canceled out.
   */
  if( wg/taud >= 20.0)
  {
    x1 = -(DT/wg - wg/taud)/sqrt(2.0);
    x2 =  (DT/wg + wg/taud)/sqrt(2.0);

    St = exp( - DT*DT/wg/wg/2.0 ) / sqrt(PI) * ( 1.0/x1 * ( 1.0 - 1.0/(2*x1*x1) +  3.0/pow(2*x1*x1, 2) - 15.0/pow(2*x1*x1, 3) ) 
                                                +1.0/x2 * ( 1.0 - 1.0/(2*x2*x2) +  3.0/pow(2*x2*x2, 2) - 15.0/pow(2*x2*x2, 3)) ) ;
  }
  else
  {
    St = exp(wg*wg/2.0/taud/taud) * ( exp(-DT/taud) * erfc( -(DT/wg - wg/taud)/sqrt(2.0) ) 
                                     +exp( DT/taud) * erfc(  (DT/wg + wg/taud)/sqrt(2.0) ));
  }

  St *= 1.0/2.0 * fg * sigma*sigma;
  return St;
}


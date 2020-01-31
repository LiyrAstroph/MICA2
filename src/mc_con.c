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

#include "dnest_con.h"
#include "allvars.h"

#include "proto.h"

void *best_model_con;   /*!< best model */
void *best_model_std_con;  /*!< standard deviation of the best model */

void mc_con()
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
  strcat(argv[argc++], "/data/restart_dnest.txt");

  /* level-dependent sampling */
  {
    strcpy(argv[argc++], "-l");
  }

  if(flag_postprc == 1)
    strcpy(argv[argc++], "-p");

  mc_con_init();

  dnest_con(argc, argv);
  
  postprocess_con();
  
  if(thistask == roottask)
  {
    FILE *fp;
    char fname[200];
    int j, idx;
    double sigma, tau, alpha, syserr, tspan;
    double *pm = best_model_con;
    double *tcon_data, *fcon_data, *fecon_data, *tcon, *fcon, *fecon;
    int ncon, ncon_data;

    sprintf(fname, "%s/%s", parset.file_dir, "data/pcon.txt");
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
 
    tcon = malloc(ncon_max*5 * sizeof(double));
    fcon = malloc(ncon_max*5 * sizeof(double));
    fecon = malloc(ncon_max*5 * sizeof(double));

    for(i=0; i<nset; i++)
    {

      idx = idx_con_pm[i];
      syserr = (exp(pm[idx])-1.0)*dataset[i].con.error_mean;
      tau = exp(pm[idx+2]);
      sigma = exp(pm[idx+1]) * sqrt(tau);
      alpha = 1.0;
      
      ncon_data = dataset[i].con.n;
      tcon_data = dataset[i].con.t;
      fcon_data = dataset[i].con.f;
      fecon_data = dataset[i].con.fe;

      ncon = ncon_data * 5; 

      tspan = tcon_data[ncon_data-1] - tcon_data[0];
      for(j=0; j<ncon; j++)
        tcon[j] = (tspan + 0.1*tspan)/(ncon-1.0) * j + tcon_data[0]-0.05*tspan;

      recostruct_con_from_varmodel(sigma, tau, alpha, syserr, ncon_data, tcon_data, fcon_data, fecon_data, ncon, tcon, fcon, fecon);

      for(j=0; j<ncon; j++)
      {
        fprintf(fp, "%f %f %f\n", tcon[j], fcon[j] * dataset[i].con.scale, fecon[j] * dataset[i].con.scale);
      }
      fprintf(fp, "\n");
    }
    
    fclose(fp);

    free(tcon);
    free(fcon);
    free(fecon);
  }
  
  mc_con_end();

  /* clear up argv */
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void recostruct_con_from_varmodel(double sigma, double tau, double alpha, double syserr, 
  int ncon_data, double *tcon_data, double *fcon_data, double *fecon_data, int ncon, double *tcon, double *fcon, double *fecon)
{
  double *Larr, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq;
  int i, j, info, *ipiv;
  double *PEmat1, *PEmat2, *PEmat3, *PEmat4;

  Larr = workspace;
  ybuf = Larr + ncon_data * nq;
  y = ybuf + ncon_data;
  Cq = y + ncon_data;
  yq = Cq + nq*nq;
  yuq = yq + nq; 
  Larr_rec = yuq + ncon;

  ipiv = workspace_ipiv;

  PEmat1 = malloc(ncon * ncon_data * sizeof(double));
  PEmat2 = malloc(ncon * ncon * sizeof(double));
  PEmat3 = malloc(ncon * nq * sizeof(double));
  PEmat4 = malloc(ncon * ncon * sizeof(double));

  set_covar_Pmat_data(sigma, tau, alpha, syserr, ncon_data, tcon_data, fcon_data, fecon_data);

  set_covar_Umat(sigma, tau, alpha, ncon_data, tcon_data, ncon, tcon);

  inverse_mat(PCmat, ncon_data, &info, ipiv);

  for(i=0;i<ncon_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(tcon_data[i], j);
  }

  /* L^T x C^-1 x L */
  multiply_mat_MN(PCmat, Larr, ybuf, ncon_data, nq, ncon_data);
  multiply_mat_MN_transposeA(Larr, ybuf, Cq, nq, nq, ncon_data);

  /* L^T x C^-1 x y */
  multiply_matvec(PCmat, fcon_data, ncon_data, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yuq, nq, 1, ncon_data);

  /* (L^T x C^-1 x L)^-1 x  L^T x C^-1 x y */
  inverse_mat(Cq, nq, &info, ipiv);
  multiply_mat_MN(Cq, yuq, yq, nq, 1, nq);

  /*  L x q */
  multiply_matvec_MN(Larr, ncon_data, nq, yq, ybuf);
  for(i=0; i<ncon_data; i++)
  {
    y[i] = fcon_data[i] - ybuf[i];
  }

  /* S x C^-1 x (y - Lq) */
  multiply_matvec(PCmat, y, ncon_data, ybuf);
  multiply_matvec_MN(USmat, ncon, ncon_data, ybuf, fcon);
  
  for(i=0;i<ncon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(tcon[i], j);
  }
  multiply_matvec_MN(Larr_rec, ncon, nq, yq, yuq);

  for(i=0; i<ncon; i++)
    fcon[i] += yuq[i];
  
  /* S x C^-1 x S */
  multiply_mat_MN(USmat, PCmat, PEmat1, ncon, ncon_data, ncon_data);
  multiply_mat_MN_transposeB(PEmat1, USmat, PEmat2, ncon, ncon, ncon_data);

  /* S x C^-1 x L */
  multiply_mat_MN(PEmat1, Larr, PEmat3, ncon, nq, ncon_data);
  /* S x C^-1 x L - L */
  for(i=0; i<ncon*nq; i++)PEmat3[i] -= Larr_rec[i];
  /* (S x C^-1 x L - L) x Cq */
  multiply_mat_MN(PEmat3, Cq, PEmat1, ncon, nq, nq);
  /* (S x C^-1 x L - L) x Cq x (S x C^-1 x L - L)^T */
  multiply_mat_MN_transposeB(PEmat1, PEmat3, PEmat4, ncon, ncon, nq);

  for(i=0; i<ncon; i++)
  {
    fecon[i] = sigma * sqrt(1.0 - PEmat2[i*ncon+i] + PEmat4[i*ncon + i]);
  }

  free(PEmat1);
  free(PEmat2);
  free(PEmat3);
  free(PEmat4);
  return;
}

/*!
 *  this function does postprocess for continuum. 
 */
void postprocess_con()
{
  char posterior_sample_file[MICA_MAX_STR_LENGTH];
  double *pm, *pmstd;
  int num_ps, i, j;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);
  
  best_model_con = malloc(size_of_modeltype);
  best_model_std_con = malloc(size_of_modeltype);
  
  if(thistask == roottask)
  {
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
    pm = (double *)best_model_con;
    pmstd = (double *)best_model_std_con;
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

    for(j = 0; j<num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_con + j), 
                                             *((double *)best_model_std_con+j) ); 

    free(post_model);
    free(posterior_sample);

    memcpy(var_param, best_model_con, num_params_var*sizeof(double));
    memcpy(var_param_std, best_model_std_con, num_params_var*sizeof(double));
  }

  MPI_Bcast(var_param, num_params_var, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(var_param_std, num_params_var, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  return;
}

double prob_con_variability(const void *model)
{
  double prob = 0.0;
  int i, k, info, sign, *ipiv;
  double *pm = (double *)model;
  double tau, sigma, alpha, lndet, lndet_ICq, syserr;
  double *Larr, *ybuf, *y, *yq, *Cq, *ICq;
  double *tcon, *fcon, *fecon;
  int ncon, idx;

  Larr = workspace;
  ybuf = Larr + ncon_max;
  y = ybuf + ncon_max;
  yq = y + ncon_max;
  Cq = yq + nq;
  ICq = Cq + nq * nq;

  ipiv = workspace_ipiv;

  for(i=0;i<ncon_max;i++)
  {
    Larr[i*nq]=1.0;
  }

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    ncon = dataset[k].con.n;
    tcon = dataset[k].con.t;
    fcon = dataset[k].con.f;
    fecon = dataset[k].con.fe;

    idx = idx_con_pm[k];
    syserr = (exp(pm[idx])-1.0)*dataset[k].con.error_mean;
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]) * sqrt(tau);
    alpha = 1.0;
    
    set_covar_Pmat_data(sigma, tau, alpha, syserr, ncon, tcon, fcon, fecon);
    memcpy(IPCmat, PCmat, ncon*ncon*sizeof(double));

    inverse_mat(IPCmat, ncon, &info, ipiv); /* calculate C^-1 */

    /* calculate L^T*C^-1*L */
    multiply_mat_MN(IPCmat, Larr, ybuf, ncon, nq, ncon);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nq, nq, ncon);
    memcpy(ICq, Cq, nq*nq*sizeof(double));
    

    /* calculate L^T*C^-1*y */
    multiply_matvec(IPCmat, fcon, ncon, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, ncon);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_mat(Cq, nq, &info, ipiv);
    multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);
  
    multiply_matvec_MN(Larr, ncon, nq, ybuf, y);
    for(i=0; i<ncon; i++)
    {
      y[i] = fcon[i] - y[i];
    }
  
    /* y^T x C^-1 y */
    multiply_matvec(IPCmat, y, ncon, ybuf);
    prob += -0.5 * cblas_ddot(ncon, y, 1, ybuf, 1)/(sigma * sigma);

    lndet = lndet_mat3(PCmat, ncon, &info, &sign, ipiv) + 2.0*ncon*log(sigma);
    if(info!=0|| sign==-1)
    {
      prob = -DBL_MAX;
      printf("lndet_C %f %d!\n", lndet, sign);
      return prob;
    }
    lndet_ICq = lndet_mat3(ICq, nq, &info, &sign, ipiv) - 2.0*nq*log(sigma);
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

double prob_con_variability_semiseparable(const void *model)
{
  double prob = 0.0;
  int i, k, info, sign, *ipiv;
  double *pm = (double *)model;
  double tau, sigma, sigma2, alpha, lndet, lndet_ICq, syserr;
  double *Larr, *ybuf, *y, *yq, *Cq;
  double *tcon, *fcon, *fecon;
  int ncon, idx;
  double *W, *D, *phi;

  ipiv = workspace_ipiv;

  Larr = workspace;
  ybuf = Larr + ncon_max;
  y = ybuf + ncon_max;
  yq = y + ncon_max;
  Cq = yq + nq;
  W = Cq + nq*nq;
  D = W + ncon_max;
  phi = D + ncon_max;

  for(i=0;i<ncon_max;i++)
  {
    Larr[i*nq]=1.0;
  }

  /* iterate over all datasets */
  for(k=0; k<nset; k++)
  {
    ncon = dataset[k].con.n;
    tcon = dataset[k].con.t;
    fcon = dataset[k].con.f;
    fecon = dataset[k].con.fe;

    idx = idx_con_pm[k];
    syserr = (exp(pm[idx])-1.0)*dataset[k].con.error_mean;
    tau = exp(pm[idx+2]);
    sigma = exp(pm[idx+1]) * sqrt(tau);
    alpha = 1.0;
    
    sigma2 = sigma*sigma;
    compute_semiseparable_drw(tcon, ncon, sigma2, 1.0/tau, fecon, syserr, W, D, phi);
    
    lndet = 0.0;
    for(i=0; i<ncon; i++)
      lndet += log(D[i]);

    // Cq^-1 = L^TxC^-1xL
    multiply_mat_semiseparable_drw(Larr, W, D, phi, ncon, nq, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, Cq, nq, nq, ncon);
    
    // L^TxC^-1xy
    multiply_matvec_semiseparable_drw(fcon, W, D, phi, ncon, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, ncon);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_symat_lndet(Cq, nq, &lndet_ICq, &info, ipiv);
    multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);
  
    multiply_matvec_MN(Larr, ncon, nq, ybuf, y);
    for(i=0; i<ncon; i++)
    {
      y[i] = fcon[i] - y[i];
    }
  
    /* y^T x C^-1 y */
    multiply_matvec_semiseparable_drw(y, W, D, phi, ncon, sigma2, ybuf);
    prob += -0.5 * cblas_ddot(ncon, y, 1, ybuf, 1);
    prob += -0.5*lndet - 0.5*lndet_ICq;
  }
  
  return prob;
}

int mc_con_init()
{

  sprintf(dnest_options_file, "%s/src/OPTIONSCON", parset.file_dir);

  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  return 0;
}

int mc_con_end()
{

  free(best_model_con);
  free(best_model_std_con);

  return 0;
}

/*!
 * this function sets the covariance matrix at data time points 
 */
void set_covar_Pmat_data(double sigma, double tau, double alpha, double syserr, int n, double *t, double *f, double *fe)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n; i++)
  { 
    t1 = t[i];

    for(j=0; j<i; j++)
    {
      t2 = t[j];
      PSmat[i*n+j] = exp (- pow (fabs(t1-t2) / tau, alpha));
      PSmat[j*n+i] = PSmat[i*n+j];

      PNmat[i*n+j] = PNmat[j*n+i] = 0.0;

      PCmat[i*n+j] = PCmat[j*n+i] = PSmat[i*n+j];
    }

    PSmat[i*n+i] = 1.0;
    PNmat[i*n+i] = (fe[i]*fe[i] + syserr*syserr)/(sigma * sigma);
    PCmat[i*n+i] = PSmat[i*n+i] + PNmat[i*n+i];
  }
  return;
}

/*!
 * this function sets the covariance matrix at time of data points and reconstruction points
 */
void set_covar_Umat(double sigma, double tau, double alpha, int ncon_data, double *tcon_data, int ncon, double *tcon)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<ncon; i++)
  {
    t1 = tcon[i];
    for(j=0; j<ncon_data; j++)
    {
      t2 = tcon_data[j];
      USmat[i*ncon_data+j] =  exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}

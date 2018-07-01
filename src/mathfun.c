/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <lapacke.h>

#include "allvars.h"
#include "proto.h"
/*
 * note that for a matrix A(mxn) in the row-major order, the first dimension is n. 
 */

/* C(nxn) = A(nxn)*B(nxn) */
void multiply_mat(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}
/* C(nxn) = A^T(nxn)*B(nxn) */
void multiply_mat_transposeA(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}
/* C(nxn) = A(nxn)*B^T(nxn) */
void multiply_mat_transposeB(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}
/* C(m*n) = A(m*k) * B(k*n) */
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, k, b, n, 0.0f, c, n);
}
/* C(m*n) = A^T(m*k) * B(k*n) */
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, n);
}
/* C(m*n) = A(m*k) * B^T(k*n) */
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0f
                             , a, k, b, k, 0.0f, c, n);
}
/* y(n) = A(nxn)* x(n) */
void multiply_matvec(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}
/* y(n) = A^T(nxn)* x(n) */
void multiply_matvec_transposeA(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}
/* y(m) = A(m, n) * x(n) */
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}
/* A(mxm)^-1 * B(mxn), store the output in B
 * note that A will be changed on exit. */
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n)
{
  int * ipiv, info;
  ipiv=malloc(m*sizeof(int));

  info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, a, m, ipiv);
  if(info!=0)
  {
    strcpy(str_error_exit, "multiply_mat_MN_inverseA 1.\n this usually caused by improper nc.\n increase the low limit of nc");
    error_exit(9);
  }
  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', m, n, a, m, ipiv, b, n);
  if(info!=0)
  {
    strcpy(str_error_exit, "multiply_mat_MN_inverseA 2\n this usually caused by improper nc.\n increase the low limit of nc");
    error_exit(9);
  }

  free(ipiv);
  return info;
}
/* A^-1 */
void inverse_mat(double * a, int n, int *info)
{
  int * ipiv;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  *info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat");
    error_exit(9);
  }
  *info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat");
    error_exit(9);
  }
  free(ipiv);
  return;
}
/* eigen vector and eigen values */
void eigen_sym_mat(double *a, int n, double *val, int *info)
{
    char jobz='V', uplo='U';

/* store the eigenvectors  in a by rows.
 * store the eigenvalues in val in ascending order.
 */
//    dsyev_(&jobz, &uplo, &n, a, &n, val, work, &lwork, info);

/* store the eigenvectors  in a by columns.
 * store the eigenvalues in val in ascending order.
 */
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, val);
    return;
}
/* A(nxn) = x^T(n) * x(n) */
void multiply_vec2mat(double * x, double * a, int n)
{
//  cblas_dsyr(CblasRowMajor, CblasUpper, n, 1.0f, x, 1, a, n);
  int i, j;
  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      a[i*n+j] = a[j*n+i] = x[i]*x[j];
    }
}
/* determinant of matrix A 
 * note that A is changed on exit. */
double det_mat(double *a, int n, int *info)
{
  int *ipiv;
  int i;
  double det;

  ipiv = workspace_ipiv;

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  
  if(*info!=0)
  {
    strcpy(str_error_exit, "det_mat");
    error_exit(9);
  }

  det = 1.0;
  for(i=0; i<n; i++)
  {
    det *= a[i*n+i];
    if (ipiv[i] != i+1)
    {
      det = -det;
    }
  }
  return det;
}
/* natural logarithm of determinant of A
 * note that A is changed on exit. */
double lndet_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double lndet;
  ipiv=workspace_ipiv;

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "lndet_mat");
    error_exit(9);
  }
  lndet = 0.0;
  for(i=0; i<n; i++)
  {
      lndet += log(fabs(a[i*n+i]));
  }
  return lndet;
}

/* natural logarithm of determinant of A
 * if any sign of the eigen values is negtive, return sign=-1 
 * note that A is changed on exit 
 */
double lndet_mat2(double *a, int n, int *info, int *sign)
{
  int * ipiv;
  int i, *sign_all;
  double lndet;
  //ipiv=malloc(n*sizeof(int));
  //sign_all=malloc(n*sizeof(int));
  ipiv = workspace_ipiv;
  sign_all = ipiv + n;

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "lndet_mat2");
    error_exit(9);
  }

  lndet = 0.0;
  for(i=0; i<n; i++)
  {
    lndet += log(fabs(a[i*n+i]));
    //*sign *= (a[i*n+i]>=0?1:-1);
    sign_all[i] = (a[i*n+i]>=0?1:-1);

    if(ipiv[i]!=i+1)
    {
      //printf("%e %d\n", a[i*n+i], n);
      //ipiv[ipiv[i]] = ipiv[i];
      //*sign *= -1;
      sign_all[i] *= -1;
    }
  }

  for(i=0; i<n;i++)
  {
    if(sign_all[i]==-1)
    {
      *sign=-1;      // if any sign is negtive, return sign=-1
      break;
    }
  }
  return lndet;
}
/* natural logarithm of determinant of A
 * sign of all the eigen values
 * note that A is changed on exit 
 */
double lndet_mat3(double *a, int n, int *info, int *sign)
{
  int *ipiv;
  int i;
  double lndet;
  //ipiv=malloc(n*sizeof(int));
  //sign_all=malloc(n*sizeof(int));
  ipiv = workspace_ipiv;

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "lndet_mat3");
    error_exit(9);
  }

  lndet = 0.0;
  *sign = 1;
  for(i=0; i<n; i++)
  {
    lndet += log(fabs(a[i*n+i]));
    *sign *= (a[i*n+i]>=0?1:-1);

    if(ipiv[i]!=i+1)
    {
      //printf("%e %d\n", a[i*n+i], n);
      //ipiv[ipiv[i]] = ipiv[i];
      *sign *= -1;
    }
  }

  return lndet;
}
/* Cholesky decomposition of A, i.e., A = M^T* M 
 * store M in A */
void Chol_decomp_U(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'L';
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    strcpy(str_error_exit, "Chol_decomp_U");
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
    error_exit(9);
  }
  else if (*info>0)
  {
    strcpy(str_error_exit, "Chol_decomp_U");
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
    error_exit(9);
  }
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      a[i*n+j] = 0.0;
  return;
}
void display_mat(double *a, int m, int n)
{
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%f\t", a[i*n+j]);
        }
        printf("\n");
    }
}

double ** matrix_malloc(int n1, int n2)
{
  double **mat;
  int i;

  if(!(mat = malloc(n1*sizeof(double*))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    error_exit(7);
  }

  for(i=0; i<n1; i++)
  {
    if(!(mat[i] = malloc(n2*sizeof(double))))
    {
      strcpy(str_error_exit, "matrix");
      error_exit(7);
    }
  }
  return mat;
}

double * array_malloc(int n)
{
  double *array;

  if(!(array = malloc(n*sizeof(double))))
  {
    strcpy(str_error_exit, "array");
    error_exit(7);
  }

  return array;
}

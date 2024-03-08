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
/* C(m*n) = - A(m*k) * B(k*n) */
void multiply_mat_MN_alpha(double * a, double *b, double *c, int m, int n, int k, double alpha)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha
                             , a, k, b, n, 0.0f, c, n);
}
/* C(m*n) = A^T(m*k) * B(k*n) */
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, n);
}
/* C(m*n) = A^T(m*k) * B(k*n) */
void multiply_mat_MN_transposeA_alpha_beta(double * a, double *b, double *c, int m, int n, int k
                                          ,double alpha, double beta)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha
                             , a, m, b, n, beta, c, n);
}
/* C(m*n) = A(m*k) * B^T(k*n) */
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0f
                             , a, k, b, k, 0.0f, c, n);
}
/* C(m*n) = alpha * A(m*k) * B^T(k*n) + beta * C(m*n) */
void multiply_mat_MN_transposeB_alpha_beta(double * a, double *b, double *c, int m, int n, int k, 
                                           double alpha, double beta)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha
                             , a, k, b, k, beta, c, n);
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
/* y(m) = alpha * A(m, n) * x(n) + beta * y(m) */
void multiply_matvec_MN_alpha_beta(double * a, int m, int n, double *x, double *y, double alpha, double beta)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, a, n, x, 1, beta, y, 1);
}
/* y(m) = A^T(m, n) * x(n) */
void multiply_matvec_MN_transposeA(double * a, int m, int n, double *x, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}
/* A(mxm)^-1 * B(mxn), store the output in B
 * note that A will be changed on exit. */
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n, int *ipiv)
{
  int info;

  info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, a, m, ipiv);
  if(info!=0)
  {
    strcpy(str_error_exit, "error in multiply_mat_MN_inverseA.");
    error_exit(9);
  }
  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', m, n, a, m, ipiv, b, n);
  if(info!=0)
  {
    strcpy(str_error_exit, "error in multiply_mat_MN_inverseA.");
    error_exit(9);
  }

  return info;
}
/*! A^-1
 *
 *  Note: ipiv, size at least max(1,min(m, n)). 
 *  Contains the pivot indices; for 1 ≤i≤ min(m, n), 
 *  row i was interchanged with row ipiv(i).
 */
void inverse_mat(double * a, int n, int *info, int *ipiv)
{
//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  *info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat1");
    error_exit(9);
  }
  *info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat2");
    error_exit(9);
  }
  return;
}

void inverse_mat_lndet(double * a, int n, double *lndet, int *info, int *sign, int *ipiv)
{
  int i;

//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  *info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat_lndet1");
    error_exit(9);
  }

  *lndet = 0.0;
  *sign = 1;
  for(i=0; i<n; i++)
  {
    *lndet += log(fabs(a[i*n+i]));
    *sign *= (a[i*n+i]>=0?1:-1);

    if(ipiv[i]!=i+1)
    {
      //printf("%e %d\n", a[i*n+i], n);
      //ipiv[ipiv[i]] = ipiv[i];
      *sign *= -1;
    }
  }

  *info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_mat_lndet2");
    error_exit(9);
  }

  return;
}
/*! A^-1, A is symmetric.
 *  
 * ipiv: array, size at least max(1, n). 
 * Contains details of the interchanges and the block structure of D. 
 * If ipiv[i-1] = k >0, then dii is a 1-by-1 block, and the i-th row 
 * and column of A was interchanged with the k-th row and column.
 * 
 * If uplo = 'U' and ipiv[i] =ipiv[i-1] = -m < 0, then D has a 2-by-2 
 * block in rows/columns i and i+1, and i-th row and column of A was 
 * interchanged with the m-th row and column.
 * 
 * If uplo = 'L' and ipiv[i] =ipiv[i-1] = -m < 0, then D has a 2-by-2 
 * block in rows/columns i and i+1, and (i+1)-th row and column of 
 * A was interchanged with the m-th row and column.
 */
void inverse_symat(double * a, int n, int *info, int *ipiv)
{
  int i, j;

  *info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_sysmat1");
    error_exit(9);
  }

  *info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_sysmat2");
    error_exit(9);
  }

  /* fill up the lower triangle */
  for(i=0; i<n; i++)
    for(j=0; j<i; j++)
      a[i*n+j] = a[j*n+i];

  return;
}
void inverse_symat_lndet(double * a, int n, double *lndet, int *info, int *ipiv)
{
  int i, j;

  *info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  
  *lndet = 0.0;
  for(i=0; i<n; i++)
  {
    *lndet += log(a[i*n+i]);
  }

  *info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);

  /* fill up the lower triangle */
  for(i=0; i<n; i++)
    for(j=0; j<i; j++)
      a[i*n+j] = a[j*n+i];
  return;
}
/* A^-1 */
void inverse_symat_lndet_sign(double * a, int n, double *lndet, int *info, int *sign, int *ipiv)
{
  int i, j;

  *info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_sysmat_lndet_sign1");
    //error_exit(9);
  }

  *lndet = 0.0;
  *sign = 1;
  for(i=0; i<n; i++)
  {
    *lndet += log(fabs(a[i*n+i]));
    *sign *= (a[i*n+i]>=0?1:-1);
  }

  *info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, a, n, ipiv);
  if(*info!=0)
  {
    strcpy(str_error_exit, "inverse_sysmat_lndet_sign2");
    //error_exit(9);
  }

  /* fill up the lower triangle */
  for(i=0; i<n; i++)
    for(j=0; j<i; j++)
      a[i*n+j] = a[j*n+i];

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
double det_mat(double *a, int n, int *info, int *ipiv)
{
  int i;
  double det;

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
double lndet_mat(double *a, int n, int *info, int *ipiv)
{
  int i;
  double lndet;

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
double lndet_mat2(double *a, int n, int *info, int *sign, int *ipiv)
{
  int i;
  double lndet;

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
  }

  for(i=0; i<n; i++)
  {
    *sign = (a[i*n+i]>=0?1:-1);
    if(ipiv[i]!=i+1)
    {
      *sign *= -1;
    }
    if(*sign==-1)
      break;
  }
  return lndet;
}
/* natural logarithm of determinant of A
 * sign of all the eigen values
 * note that A is changed on exit 
 */
double lndet_mat3(double *a, int n, int *info, int *sign, int *ipiv)
{
  int i;
  double lndet;

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
            printf("%e\t", a[i*n+j]);
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

/*
 * semiseparable matrix
 */
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr, double *W, double *D, double *phi)
{
  int i;
  double S, A;
  phi[0] = 0.0;
  for(i=1; i<n; i++)
  {
    phi[i] = exp(-c1 * (t[i] - t[i-1]));
  }

  S = 0.0;
  A = sigma[0]*sigma[0] + syserr*syserr + a1;
  D[0] = A;
  W[0] = 1.0/D[0];
  for(i=1; i<n; i++)
  {
    S = phi[i]*phi[i] * (S + D[i-1]*W[i-1]*W[i-1]);
    A = sigma[i]*sigma[i] + syserr*syserr + a1;
    D[i] = A - a1 * a1 * S;
    W[i] = 1.0/D[i] * (1.0 - a1*S);
  }
}
/*
 * z = C^-1 x y
 *
 * y is a vector
 */
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z)
{
  int i;
  double f, g;

  // forward substitution
  f = 0.0;
  z[0] = y[0];
  for(i=1; i<n;i++)
  {
    f = phi[i] * (f + W[i-1] * z[i-1]);
    z[i] = y[i] - a1*f;
  }

  //backward substitution
  g = 0.0;
  z[n-1] = z[n-1]/D[n-1];
  for(i=n-2; i>=0; i--)
  {
    g = phi[i+1] *(g + a1*z[i+1]);
    z[i] = z[i]/D[i] - W[i]*g;
  }
}
/*
 * Z = C^-1 x Y
 * 
 * Y is an (nxm) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0*m+j];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i*m+j] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

/*
 * Z = C^-1 x Y^T
 * 
 * Y is an (mxn) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0+j*n];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i+j*n] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

/*!
 * This function calculates A^-1 by partitioning, where A is symmetric.
 * 
 * A is partitioned into 
 *  
 *   P  ,  Q
 *   Q^T,  S
 * 
 * where P(n1xn1), Q(n1xn2), S(n2xn2), n=n1+n2.
 * 
 * The size of work is (n1*n2 + n1*n1)=n1*n.
 */
void inverse_symat_partition(double *P, double *S, double *Q, int n1, int n2, 
                             double *PN, double *SN, double *QN, double *work, 
                             int *ipiv)
{
  int i, j, info;
  double *pwork;

  /* P^-1 */
  memcpy(PN, P, (n1*n1)*sizeof(double));
  inverse_symat(PN, n1, &info, ipiv);

  /* P^-1 x Q */
  multiply_mat_MN(PN, Q, work, n1, n2, n1);

  /* Q^T x P^-1 x Q */
  multiply_mat_MN_transposeA(Q, work, SN, n2, n2, n1);

  /* (S - Q^T x P^-1 x Q)^-1; only upper triangle */
  for(i=0; i<n2; i++)
  {
    for(j=i; j<n2; j++)
    {
      SN[i*n2 + j] = S[i*n2 + j] - SN[i*n2 + j];
    }
  }
  inverse_symat(SN, n2, &info, ipiv);

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 */
  multiply_mat_MN(work, SN, QN, n1, n2, n2);

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 x (P^-1 x Q)^T */
  pwork = work+n1*n2;
  multiply_mat_MN_transposeB(QN, work, pwork, n1, n1, n2);

  for(i=0; i<n1*n2; i++)
  {
    QN[i] = -QN[i];
  }

  for(i=0; i<n1*n1; i++)
  {
    PN[i] = PN[i] + pwork[i];
  }

  return;
}

/*!
 * This function calculates A^-1 by partitioning and given P^-1, 
 * where A is symmetric.
 * 
 * A is partitioned into 
 *  
 *   P  ,  Q
 *   Q^T,  S
 * 
 * where P(n1xn1), Q(n1xn2), S(n2xn2), n=n1+n2.
 * 
 * The size of work is is (n1*n2 + n1*n1)=n1*n.
 */
void inverse_symat_partition_inv(double *Pinv, double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *work, int *ipiv)
{
  int i, j, info;
  double *pwork;

  /* P^-1 x Q */
  multiply_mat_MN(Pinv, Q, work, n1, n2, n1);

  /* Q^T x P^-1 x Q */
  multiply_mat_MN_transposeA(Q, work, SN, n2, n2, n1);

  /* (S - Q^T x P^-1 x Q)^-1; only upper triangle */
  for(i=0; i<n2; i++)
  {
    for(j=i; j<n2; j++)
    {
      SN[i*n2 + j] = S[i*n2 + j] - SN[i*n2 + j];
    }
  }
  inverse_symat(SN, n2, &info, ipiv);

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 */
  multiply_mat_MN(work, SN, QN, n1, n2, n2);

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 x (P^-1 x Q)^T */
  pwork = work+n1*n2;
  multiply_mat_MN_transposeB(QN, work, pwork, n1, n1, n2);

  for(i=0; i<n1*n2; i++)
  {
    QN[i] = -QN[i];
  }

  for(i=0; i<n1*n1; i++)
  {
    PN[i] = Pinv[i] + pwork[i];
  }

  return;
}

/*
 * calculate matrix inverse by partition
 * ref: Numerical Recipes, Chapter 2, P70
 */
void inverse_symat_lndet_partition_inv(double *Pinv, double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *lndet, double *work, int *ipiv)
{
  int i, j, info;
  double *pwork;

  /* P^-1 x Q */
  multiply_mat_MN(Pinv, Q, work, n1, n2, n1);

  /* Q^T x P^-1 x Q */
  multiply_mat_MN_transposeA(Q, work, SN, n2, n2, n1);

  /* (S - Q^T x P^-1 x Q)^-1; only upper triangle*/
  for(i=0; i<n2; i++)
  {
    for(j=i; j<n2; j++)
    {
      SN[i*n2 + j] = S[i*n2 + j] - SN[i*n2 + j];
    }
  }
  /* note that lndet = lndet(SN^-1) */
  inverse_symat_lndet(SN, n2, lndet, &info, ipiv); 

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 */
  multiply_mat_MN(work, SN, QN, n1, n2, n2);

  /* (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 x (P^-1 x Q)^T */
  pwork = work+n1*n2;
  multiply_mat_MN_transposeB(QN, work, pwork, n1, n1, n2);

  for(i=0; i<n1*n2; i++)
  {
    QN[i] = -QN[i];
  }

  for(i=0; i<n1*n1; i++)
  {
    PN[i] = Pinv[i] + pwork[i];
  }

  return;
}
/*
 * calculate matrix inverse by partition
 * ref: Numerical Recipes, Chapter 2, P70
 * use the functionality such as C = a * AxB + b * C in Lapack.
 */
void inverse_symat_lndet_partition_inv_fast(double *Pinv, double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *lndet, double *work, int *ipiv)
{
  int info;

  /* P^-1 x Q */
  multiply_mat_MN(Pinv, Q, work, n1, n2, n1);

  /* S - Q^T x P^-1 x Q */
  // memcpy(SN, S, n2*n2*sizeof(double)); /* note no need to copy S, since SN=S */
  multiply_mat_MN_transposeA_alpha_beta(Q, work, SN, n2, n2, n1, -1.0, 1.0);

  /* note that lndet = lndet(SN^-1) */
  inverse_symat_lndet(SN, n2, lndet, &info, ipiv); 

  /* -(P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 */
  multiply_mat_MN_alpha(work, SN, QN, n1, n2, n2, -1.0);

  /* P^-1 + (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 x (P^-1 x Q)^T */
  memcpy(PN, Pinv, n1*n1*sizeof(double));
  multiply_mat_MN_transposeB_alpha_beta(QN, work, PN, n1, n1, n2, -1.0, 1.0);

  return;
}

/*!
 * This function calculates A^-1 by iterative partitioning,
 * where A is symmetric.
 * 
 * A is partitioned into 
 *  
 *   An  ,  Qn
 *   Qn^T,  Sn
 * 
 * and An is partitioned into 
 * 
 *   An-1  ,  Qn-1
 *   Qn-1^T,  Sn-1
 * 
 * and so on ... in total n iterations.
 * 
 * narr stores the size of A1, S1， S2，..., Sn. note the size 
 * of narr is n+1.
 * 
 * The size of work_inv is about 3*nt*nt.
 * 
 * On exit, calculate lndet(A).
 */
void inverse_symat_partition_iter(double *A, int nt, int *narr, int np, double *lndet, 
                                  double *work_inv, int *ipiv)
{
  int i, j, k, ni, nq, info;
  double *Ai, *Qi, *Si, *ANi, *QNi, *SNi, *pwork, lndet_SN;

  ni = narr[0];
  Ai = work_inv;

  /* A1; only upper triangle */
  for(i=0; i<ni; i++)
  {
    for(j=i; j<ni; j++)
    {
      Ai[i*ni + j] = A[i*nt + j];
    }
  }
  inverse_symat_lndet(Ai, ni, lndet, &info, ipiv);
  
  for(k=1; k<=np; k++)
  {
    nq = narr[k];
    
    Qi = Ai + nt*nt;
    Si = Qi + ni*nq;
  
    ANi = Si + nq*nq;
    QNi = ANi + ni*ni;
    SNi = QNi + ni*nq;

    pwork = SNi + nq*nq;
  
    for(i=0; i<ni; i++)
    {
      for(j=0; j<nq; j++)
      {
        Qi[i*nq+j] = A[i*nt + (ni+j)];
      }
    }
  
    for(i=0; i<nq; i++)
    {
      for(j=0; j<nq; j++)
      {
        Si[i*nq+j] = A[(ni+i)*nt+(ni+j)];
      }
    }
  
    inverse_symat_lndet_partition_inv_fast(Ai, Si, Qi, ni, nq, ANi, SNi, QNi, &lndet_SN, pwork, ipiv);
    *lndet += lndet_SN; /* lndet_SN = -lndet(SN) */
    
    /* last run, directly save to A */
    if(k==np) Ai = A;
    
    /* new Ai */
    for(i=0; i<ni; i++)
    {
      /* diagonal */
      for(j=0; j<ni; j++)
      {
        Ai[i*(ni+nq) + j] = ANi[i*ni + j];
      }

      /* upper right and lower left corner */
      for(j=0; j<nq; j++)
      {
        Ai[i*(ni+nq) + (ni+j)] = Ai[(ni+j)*(ni+nq) + i] = QNi[i*nq+j];
      }
    }

    for(i=0; i<nq; i++)
    {
      for(j=0; j<nq; j++)
      {
        Ai[(ni+i)*(ni+nq) + (ni+j)] = SNi[i*nq + j];
      }
    }
    ni += nq;
  }
  
  //memcpy(A, Ai, ni*ni*sizeof(double));
  return;
}

void inverse_symat_lndet_partition_inv_semiseparable(double *Pinv, double *W, double *D, 
                                 double *phi, double a1,  
                                 double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *lndet, double *work, int *ipiv)
{
  int i, j, info;

  /* P^-1 x Q */
  multiply_mat_semiseparable_drw(Q, W, D, phi, n1, n2, a1, work);
  //multiply_mat_MN(Pinv, Q, work, n1, n2, n1);

  /* S - Q^T x P^-1 x Q */
  // memcpy(SN, S, n2*n2*sizeof(double));  /* note no need to copy S, since SN=S */
  multiply_mat_MN_transposeA_alpha_beta(Q, work, SN, n2, n2, n1, -1.0, 1.0);

  /* note that lndet = lndet(SN^-1) */
  inverse_symat_lndet(SN, n2, lndet, &info, ipiv); 

  /* -(P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 */
  multiply_mat_MN_alpha(work, SN, QN, n1, n2, n2, -1.0);

  /* P^-1 + (P^-1 x Q) x (S - Q^T x P^-1 x Q)^-1 x (P^-1 x Q)^T */
  memcpy(PN, Pinv, n1*n1*sizeof(double));
  multiply_mat_MN_transposeB_alpha_beta(QN, work, PN, n1, n1, n2, -1.0, 1.0);

  return;
}

/*!
 *  calculate A^-1 and lndet(A).
 * 
 *  A = LxDxL^T, L = I + tril(UxW^T), D is a diagonal matrix.
 * 
 *  M = LxD^1/2,  A = MxM^T,  A^-1 = (M^T)^-1xM^-1.
 */
void inverse_semiseparable(double *t, int n, double a1, double c1, double *sigma, 
                           double syserr, double *W, double *D, double *phi,
                           double *A, double *lndet)
{
  int i, j;

  compute_semiseparable_drw(t, n, a1, c1, sigma, syserr, W, D, phi);
  
  *lndet = 0.0;
  for(i=0; i<n; i++)
  {
    A[i*n + i] = 1.0 * sqrt(D[i]);
    for(j=0; j<i; j++)
    {
      A[i*n + j] = a1 * (exp(-c1*(t[i]-t[j]))*W[j]) * sqrt(D[j]);
    }

    *lndet += log(D[i]);
  }
  /* inverse of a real upper or lower triangular matrix */
  LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'N', n, A, n);
  
  /* L^T x L */
  LAPACKE_dlauum(LAPACK_ROW_MAJOR, 'L', n, A, n);

  /* fill up upper triangle */
  for(i=0; i<n; i++)
    for(j=i+1; j<n; j++)
      A[i*n+j] = A[j*n+i];
  return;
}

/*!
 *  calculate A^-1 and lndet(A) by iterative partitioning.
 * 
 *  A = LxDxL^T, L = I + tril(UxW^T), D is a diagonal matrix.
 * 
 *  M = LxD^1/2,  A = MxM^T,  A^-1 = (M^T)^-1xM^-1.
 */
void inverse_semiseparable_iter(double *t, int n, double a1, double c1, double *sigma, 
                           double syserr, double *W, double *D, double *phi,
                           double *A, int nt, int *narr, int np, double *lndet,
                           double *work_inv, int *ipiv)
{
  int i, j, k, ni, nq;
  double *Ai, *ANi, *QNi, *SNi, *pwork, lndet_SN;
 
  ni = narr[0];
  Ai = work_inv;

  inverse_semiseparable(t, ni, a1, c1, sigma, syserr, W, D, phi, Ai, lndet);
  
  for(k=1; k<=np; k++)
  {
    nq = narr[k];
    
    ANi = Ai + nt*nt;
    QNi = ANi + ni*ni;
    SNi = QNi + ni*nq;

    pwork = SNi + nq*nq;
  
    for(i=0; i<ni; i++)
    {
      for(j=0; j<nq; j++)
      {
        QNi[i*nq+j] = A[i*nt + (ni+j)];
      }
    }
  
    for(i=0; i<nq; i++)
    {
      for(j=0; j<nq; j++)
      {
        SNi[i*nq+j] = A[(ni+i)*nt+(ni+j)];
      }
    }
    
    /* for k=1, use fast calculations of semiseparable matrice */
    if(k==1)
    {
      inverse_symat_lndet_partition_inv_semiseparable(Ai, W, D, phi, a1, SNi, QNi, ni, nq, 
                                                      ANi, SNi, QNi, &lndet_SN, pwork, ipiv);
    }
    else
    {
      inverse_symat_lndet_partition_inv_fast(Ai, SNi, QNi, ni, nq, 
                                             ANi, SNi, QNi, &lndet_SN, pwork, ipiv);
    }
    *lndet += lndet_SN; /* lndet_SN = -lndet(SN) */
    
    /* last run, directly save to A */
    if(k==np) Ai = A;

    /* new Ai */
    for(i=0; i<ni; i++)
    {
      /* diagonal */
      for(j=0; j<ni; j++)
      {
        Ai[i*(ni+nq) + j] = ANi[i*ni + j];
      }

      /* upper right and lower left corner */
      for(j=0; j<nq; j++)
      {
        Ai[i*(ni+nq) + (ni+j)] = Ai[(ni+j)*(ni+nq) + i] = QNi[i*nq+j];
      }
    }

    for(i=0; i<nq; i++)
    {
      for(j=0; j<nq; j++)
      {
        Ai[(ni+i)*(ni+nq) + (ni+j)] = SNi[i*nq + j];
      }
    }
    ni += nq;
  }
  
  //memcpy(A, Ai, ni*ni*sizeof(double));
  return;
}
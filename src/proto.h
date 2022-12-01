/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef MICA_PROTO_H
#define MICA_PROTO_H

#include "allvars.h"

void begin_run();
void end_run();
void set_argv(int flag_pp, int flag_de);

/* mpi */
void set_mpi();

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);


void set_param_file(char *fname);
void set_parset(PARSET *ps);
void free_memory();
void allocate_memory();
int read_parset();
int read_data();
void scale_con_line();
void cal_mean_error();
void get_num_particles(char *fname);

void init();

void print_version();

/* for con */
double dnest_con(int argc, char **argv);
void set_par_range_con();
void mc_con();
int mc_con_init();
int mc_con_end();
double prob_con_variability(const void *model);
double prob_con_variability_semiseparable(const void *model);
void set_covar_Pmat_data(double sigma, double tau, double alpha, double syserr, int n, double *t, double *f, double *fe);
void set_covar_Umat(double sigma, double tau, double alpha, int ncon_data, double *tcon_data, int ncon, double *tcon);
void postprocess_con();
void output_reconstrction();
void output_reconstrction2();
void recostruct_con_from_varmodel(double sigma, double tau, double alpha, double syserr, 
  int ncon_data, double *tcon_data, double *fcon_data, double *fecon_data, int ncon, double *tcon, double *fcon, double *fecon);
int write_options_con(char *fname);


/* for line */
double dnest_line(int argc, char **argv);
void set_par_range_line();
void print_para_names_line();
void set_idx_line_pm();
void mc_line();
int mc_line_init();
int mc_line_end();
double prob_line_variability(const void *model);
double prob_line_variability2(const void *model);
double prob_line_variability3(const void *model);
double prob_line_variability4(const void *model);
void postprocess_line();

double (*Sll)(double t1, double t2, const void *model, int nds, int nls);
void (*Sll_array)(double *tline, int nline, const void *model, int nds, int nls, double *Smat);
double (*Sll2)(double t1, double t2, const void *model, int nds, int nls1, int nls2);
void (*Sll2_array)(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat);
double (*Slc)(double tcon, double tline, const void *model, int nds, int nls);
void (*Slc_array)(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat);

double Sll_gauss(double t1, double t2, const void *model, int nds, int nls);
void Sll_array_gauss(double *tline, int nline, const void *model, int nds, int nls, double *Smat);
double Sll2_gauss(double t1, double t2, const void *model, int nds, int nls1, int nls2);
void Sll2_array_gauss(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat);
double Slc_gauss(double tcon, double tline, const void *model, int nds, int nls);
void Slc_array_gauss(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat);

double Sll_tophat(double t1, double t2, const void *model, int nds, int nls);
void Sll_array_tophat(double *tline, int nline, const void *model, int nds, int nls, double *Smat);
double Sll2_tophat(double t1, double t2, const void *model, int nds, int nls1, int nls2);
void Sll2_array_tophat(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat);
double Slc_tophat(double tcon, double tline, const void *model, int nds, int nls);
void Slc_array_tophat(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat);

void set_covar_Pmat_data_line(const void *model, int k);
void set_covar_Pmat_data_line_array(const void *model, int k);
void set_covar_Umat_line(const void *model, int k, int *nall, double *tall);
void set_covar_Amat_line(const void *model, int k, int *nall, double *tall);
void recostruct_line_from_varmodel(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall);
void recostruct_line_from_varmodel2(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall);
int check_gauss_center(int which, int *igau);
inline int check_gauss_positivity(int which);
int write_options_line(char *fname);

// error.c
void error_exit(int);
int mica_cmp(const void * a, const void * b);

//line decomponent
void decompose_single_component(const void *model, int nds, int *nall, double *tall, double *fall, double *feall, double *yqall, int kgau);
double (*Slc_single)(double tcon, double tline, const void *model, int nds, int nls, int kgau);
double (*Sll2_single)(double t1, double t2, const void *model, int nds, int nls1, int nls2, int kgau);
double (*Sll_single)(double t1, double t2, const void *model, int nds, int nls, int kgau);
double (*Sll_single2)(double t1, double t2, const void *model, int nds, int nls, int kgau);
void set_covar_Umat_line_single(const void *model, int nds, int *nall, double *tall, int kgau);
void set_covar_Amat_line_single(const void *model, int nds, int *nall, double *tall, int kgau);
void output_decompose_line();

double Slc_single_gauss(double tcon, double tline, const void *model, int nds, int nls, int kgau);
double Sll2_single_gauss(double t1, double t2, const void *model, int nds, int nls1, int nls2, int kgau);
double Sll_single_gauss(double t1, double t2, const void *model, int nds, int nls, int kgau);
double Sll_single2_gauss(double t1, double t2, const void *model, int nds, int nls, int kgau);

double Slc_single_tophat(double tcon, double tline, const void *model, int nds, int nls, int kgau);
double Sll2_single_tophat(double t1, double t2, const void *model, int nds, int nls1, int nls2, int kgau);
double Sll_single_tophat(double t1, double t2, const void *model, int nds, int nls, int kgau);
double Sll_single2_tophat(double t1, double t2, const void *model, int nds, int nls, int kgau);

// mathfunc.c
void inverse_mat(double *a, int n, int *info, int *ipiv);
void inverse_mat_lndet(double * a, int n, double *lndet, int *info, int *sign, int *ipiv);
void inverse_symat(double *a, int n, int *info, int *ipiv);
void inverse_symat_lndet(double *a, int n, double *lndet, int *info, int *ipiv);
void inverse_symat_lndet_sign(double *a, int n, double *lndet, int *info, int *sign, int *ipiv);
double det_mat(double *a, int n, int *info, int *ipiv);
double lndet_mat(double *a, int n, int *info, int *ipiv);
double lndet_mat2(double *a, int n, int *info, int *sign, int *ipiv);
double lndet_mat3(double *a, int n, int *info, int *sign, int *ipiv);
void display_mat(double *a, int m, int n);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k);
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n, int *ipiv);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr,  double *W, double *D, double *phi);
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z);
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);
void inverse_symat_partition(double *P, double *S, double *Q, int n1, int n2, 
                             double *PN, double *SN, double *QN, double *work, int *ipiv);

void inverse_symat_partition_inv(double *Pinv, double *S, double *Q, int n1, int n2, 
                             double *PN, double *SN, double *QN, double *work, int *ipiv);


void inverse_symat_lndet_partition_inv(double *Pinv, double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *lndet, double *work, int *ipiv);

void inverse_symat_lndet_partition_inv_semiseparable(double *Pinv, double *W, double *D, 
                                 double *phi, double a1,  
                                 double *S, double *Q, int n1, 
                                 int n2, double *PN, double *SN, double *QN, 
                                 double *lndet, double *work, int *ipiv);

void inverse_symat_partition_iter(double *A, int nt, int *narr, int np, double *lndet,
                             double *work_inv, int *ipiv);

void inverse_semiseparable(double *t, int n, double a1, double c1, double *sigma, 
                           double syserr, double *W, double *D, double *phi,
                           double *A, double *lndet);


void inverse_semiseparable_iter(double *t, int n, double a1, double c1, double *sigma, 
                           double syserr, double *W, double *D, double *phi,
                           double *A, int nt, int *narr, int np, double *lndet,
                           double *work_inv, int *ipiv);

#endif
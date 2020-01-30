/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

void begin_run();
void end_run();

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);


void free_memory();
void allocate_memory();
int read_parset();
int read_data();
void scale_con_line();
void cal_mean_error();
void get_num_particles(char *fname);
void get_posterior_sample_file(char *fname, char *samplefile);

void init();


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
void output_reconstrction(const void *model);
void recostruct_con_from_varmodel(double sigma, double tau, double alpha, double syserr, 
  int ncon_data, double *tcon_data, double *fcon_data, double *fecon_data, int ncon, double *tcon, double *fcon, double *fecon);


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
void postprocess_line();
double Sll(double t1, double t2, const void *model, int nds, int nls);
void Sll_array(double *tline, int nline, const void *model, int nds, int nls, double *Smat);
double Sll2(double t1, double t2, const void *model, int nds, int nls1, int nls2);
void Sll2_array(double *tline1, int nline1, double *tline2, int nline2, const void *model, int nds, int nls1, int nls2, double *Smat);
double Slc(double tcon, double tline, const void *model, int nds, int nls);
void Slc_array(double *tcon, int ncon, double *tline, int nline, const void *model, int nds, int nls, double *Smat);
void set_covar_Pmat_data_line(const void *model, int k);
void set_covar_Pmat_data_line_array(const void *model, int k);
void set_covar_Umat_line(const void *model, int k, int *nall, double *tall);
void set_covar_Amat_line(const void *model, int k, int *nall, double *tall);
void recostruct_line_from_varmodel(const void *model, int nds, int *nall, double *tall, double *fall, double *feall);
void recostruct_line_from_varmodel2(const void *model, int nds, int *nall, double *tall, double *fall, double *feall);
int check_gauss_center(int which, int *igau);

// error.c
void error_exit(int);
int mica_cmp(const void * a, const void * b);

// mathfunc.c
void inverse_mat(double *a, int n, int *info, int *ipiv);
void inverse_mat_lndet(double * a, int n, double *lndet, int *info, int *sign);
void inverse_symat(double *a, int n, int *info);
void inverse_symat_lndet(double *a, int n, double *lndet, int *info, int *sign);
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
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n);
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

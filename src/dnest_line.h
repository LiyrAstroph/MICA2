/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_LINE_H

#include <stdbool.h>

/* size of model type, defined in dnest */
extern int size_of_modeltype;

/* functions */
void from_prior_line(void *model);
void print_particle_line(FILE *fp, const void *model);
double log_likelihoods_cal_line(const void *model);
double log_likelihoods_cal_initial_line(const void *model);
double log_likelihoods_cal_restart_line(const void *model);
double perturb_line(void *model);
int get_num_params_line();
void restart_action_line(int iflag);

double log_likelihoods_cal_line_exam(const void *model);


void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
int (*get_num_params)();
void (*restart_action)(int iflag);
#endif

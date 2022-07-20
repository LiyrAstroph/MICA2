/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_LINE_H

#define _DNEST_LINE_H

#include <stdbool.h>

/* functions */
void from_prior_line(void *model);
void print_particle_line(FILE *fp, const void *model);
double log_likelihoods_cal_line(const void *model);
double log_likelihoods_cal_initial_line(const void *model);
double log_likelihoods_cal_restart_line(const void *model);
double perturb_line_prior0(void *model);
double perturb_line_prior1(void *model);
void restart_action_line(int iflag);
double log_likelihoods_cal_line_exam(const void *model);
void read_particle_line(FILE *fp, void *model);
void kill_action_line(int i, int i_copy);
void accept_action_line();

#endif

/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_CON_H

#define _DNEST_CON_H

#include <stdbool.h>
#include <stdio.h>

/* functions */
void from_prior_con(void *model);
void print_particle_con(FILE *fp, const void *model);
double log_likelihoods_cal_con(const void *model);
double log_likelihoods_cal_initial_con(const void *model);
double log_likelihoods_cal_restart_con(const void *model);
double perturb_con(void *model);
void restart_action_con(int iflag);
double log_likelihoods_cal_con_exam(const void *model);
void read_particle_con(FILE *fp, void *model);
void kill_action_con(int i, int i_copy);
void accept_action_con();

#endif

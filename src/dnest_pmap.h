/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_PMAP_H

#define _DNEST_PAMP_H

#include <stdbool.h>
#include <stdio.h>

/* functions */
void from_prior_pmap(void *model);
void print_particle_pmap(FILE *fp, const void *model);
double log_likelihoods_cal_pmap(const void *model);
double log_likelihoods_cal_initial_pmap(const void *model);
double log_likelihoods_cal_restart_pmap(const void *model);
double perturb_pmap(void *model);
void restart_action_pmap(int iflag);
double log_likelihoods_cal_pmap_exam(const void *model);
void read_particle_pmap(FILE *fp, void *model);
void kill_action_pmap(int i, int i_copy);
void accept_action_pmap();

#endif

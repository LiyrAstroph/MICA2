/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_VMAP_H

#define _DNEST_VMAP_H

#include <stdbool.h>
#include <stdio.h>

/* functions */
void from_prior_vmap(void *model);
void print_particle_vmap(FILE *fp, const void *model);
double log_likelihoods_cal_vmap(const void *model);
double log_likelihoods_cal_initial_vmap(const void *model);
double log_likelihoods_cal_restart_vmap(const void *model);
double perturb_vmap_prior0(void *model);
double perturb_vmap_prior1(void *model);
void restart_action_vmap(int iflag);
void read_particle_vmap(FILE *fp, void *model);
void kill_action_vmap(int i, int i_copy);
void accept_action_vmap();

#endif

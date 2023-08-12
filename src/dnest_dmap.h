/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_DMAP_H

#define _DNEST_DAMP_H

#include <stdbool.h>
#include <stdio.h>

/* functions */
void from_prior_dmap(void *model);
void print_particle_dmap(FILE *fp, const void *model);
double log_likelihoods_cal_dmap(const void *model);
double log_likelihoods_cal_initial_dmap(const void *model);
double log_likelihoods_cal_restart_dmap(const void *model);
double perturb_dmap(void *model);
void restart_action_dmap(int iflag);
void read_particle_dmap(FILE *fp, void *model);
void kill_action_dmap(int i, int i_copy);
void accept_action_dmap();

#endif

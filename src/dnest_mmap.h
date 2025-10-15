/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#ifndef _DNEST_MMAP_H

#define _DNEST_MMAP_H

#include <stdbool.h>
#include <stdio.h>

/* functions */
void from_prior_mmap(void *model);
void print_particle_mmap(FILE *fp, const void *model);
double log_likelihoods_cal_mmap(const void *model);
double log_likelihoods_cal_initial_mmap(const void *model);
double log_likelihoods_cal_restart_mmap(const void *model);
double perturb_mmap_prior0(void *model);
double perturb_mmap_prior1(void *model);
void restart_action_mmap(int iflag);
double log_likelihoods_cal_mmap_exam(const void *model);
void read_particle_mmap(FILE *fp, void *model);
void kill_action_mmap(int i, int i_copy);
void accept_action_mmap();

#endif

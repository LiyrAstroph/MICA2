/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_erf.h>

#include "dnest.h"
#include "dnest_line.h"
#include "allvars.h"

#include "proto.h"

void output_decompose_line_parallel_mmap()
{
  printf("Decomposition for mmap mode is still in developing...\n");
  return;
}
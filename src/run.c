/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <string.h>

#include "proto.h"
#include "allvars.h"

#include "proto.h"

void begin_run()
{
  //strcat(parset.param_file, "src/param");

  read_parset();
  read_data();

  init();

  mc_con();

  mc_line();
}

/*!
 * This function frees the memory and ends the run.
 */
void end_run()
{
  free_memory();
}

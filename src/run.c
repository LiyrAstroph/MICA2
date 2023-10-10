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

  /* do not need to run mc_con if only postprocessing */
  if(flag_postprc!=1 && flag_decomp != 1 && parset.model != vmap)
  {
    mc_con();
  }
  
  if(parset.model == gmodel)
  {
    mc_line();
  }
  else if (parset.model == pmap)
  {
    mc_pmap();
  }
  else if (parset.model == vmap)
  {
    mc_vmap();
  }
}

/*!
 * This function frees the memory and ends the run.
 */
void end_run()
{
  free_memory();
}

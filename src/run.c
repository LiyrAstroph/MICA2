/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <string.h>
#include <sys/stat.h>
#include <stdio.h>

#include "proto.h"
#include "allvars.h"

#include "proto.h"

extern int stat(const char *path, struct stat *buf);

void begin_run()
{
  //strcat(parset.param_file, "src/param");
  
  /* check data/ and param/ directories */
  if(thistask == roottask)
  {
    struct stat st;
    if(stat("data/", &st) != 0 || !S_ISDIR(st.st_mode))
    {
      fprintf(stderr, "Error: data/ directory not found. Please create it.\n");
      exit(1);
    }
    if(stat("param/", &st) != 0 || !S_ISDIR(st.st_mode))
    {
      fprintf(stderr, "Error: param/ directory not found. Please create it.\n");
      exit(1);
    }
  }

  read_parset();
  read_data();

  init();

  /* do not need to run mc_con if only postprocessing */
  if(flag_postprc!=1 && flag_decomp != 1 && parset.model != vmap)
  {
    /* if restart, only postprocess continuum to get the parameters */
    if(flag_restart == 1)
    {
      flag_postprc = 1;
      mc_con();
      flag_postprc = 0;
    }
    else 
    {
      mc_con();
    }
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

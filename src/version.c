/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>

#include "version.h"
#include "proto.h"

void print_version()
{
  printf("\e[1;31m" "MICA Version: %d.%d.%d\n" "\e[0m", MICA_MAJOR_VERSION, MICA_MINOR_VERSION, MICA_PATCH_VERSION);
  printf("Yan-Rong Li, liyanrong@ihep.ac.cn\n");
  return;
}
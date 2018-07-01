/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

void error_exit(int code)
{
  switch(code)
  {
    case 1:
      printf("# no parameter file.\n");
      exit(1);
      break;

    case 2:
      printf("# cannot open file %s.\n", str_error_exit);
      exit(2);
      break;

    case 3:
      printf("# cannot recognize the parameter %s.\n", str_error_exit);
      exit(3);
      break;

    case 4:
      printf("# file %s does not exist.\n", str_error_exit);
      exit(4);
      break;

    case 5:
      printf("# end of parameter file %s.\n", parset.param_file);
      exit(5);
      break;

    case 6:
      printf("# error in reading file %s.\n", str_error_exit);
      exit(6);
      break;

    case 7:
      printf("# error in allocate memeory for %s.\n", str_error_exit);
      exit(7);
      break;

    case 8:
      printf("# data file %s exceeds the largest allowed size.\n", str_error_exit);
      exit(8);
      break;

    case 9:
      printf("# error in matrix calculating %s.\n", str_error_exit);
      exit(9);
      break;

    case 10:
      printf("# length of mcmc chain is smaller than nbuilt in %s.\n", str_error_exit);
      exit(10);
      break;
    case 11:
      printf("# incorrect configuration in the parameter %s.\n", str_error_exit);
      exit(11);
      break;

    case 12:
      printf("# %s.\n", str_error_exit);
      exit(12);
      break;

    default:
      printf("# cannot recognize the error code.\n");
      exit(-1);
  }
}
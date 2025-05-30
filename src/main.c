/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */


#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <mpi.h>
#include <string.h>

#include "mygetopt.h"
#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  double t0=0.0, t1=0.0, dt;
  int opt;
  extern int my_optind, my_opterr, my_optopt;
  extern char *my_optarg;
  extern int my_getopt(int argc, char *const *argv, const char *options);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);

  if(thistask == roottask)
  {
    t0 = second();
    printf("===============MICA2==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
    
    flag_postprc = 0;
    flag_end = 0;
    flag_restart = 0;
    flag_decomp = 0;
    flag_examine = 0;
    flag_load_prior = 0;
    flag_para_name = 0;
    flag_postsample = 0;
    
    my_opterr = 0; /* supress getopt error */
    my_optind = 0; /* optind=0 reset getopt */

    while( (opt = my_getopt(argc, argv, "pvrdel:ns")) != -1)
    {
      switch(opt)
      {
        case 'p':
          flag_postprc = 1;
          break;
        
        case 'v':
          flag_end = 1;
          print_version();
          break;
        
        case 'r':
          flag_restart = 1;
          break;
        
        case 'd':
          flag_postprc = 1;
          flag_decomp = 1;
          break;
        
        case 'e':
          flag_examine = 1;
          break;
        
        case 'l': /* Load parameter prior */
          flag_load_prior = 1;
          strcpy(prior_file, my_optarg);
          printf("# Load parameter prior from %s.\n", prior_file);
          break;
        
        case 'n': /* print parameter names */
          printf("# Print parameter name.\n");
          flag_para_name = 1;
          break;
        
        case 's': /* only calculate posterior sample */
          printf("# Calculate posterior sample.\n");
          flag_postsample = 1;
          break;
        
        case '?':
          printf("# Incorrect option -%c %s.\n", my_optopt, my_optarg);
          exit(0);
          break;
      }
    }

    if(my_optind < argc ) // parameter file is specified 
      strcpy(parset.param_file, argv[my_optind]); /* copy input parameter file */
    else
    {
      fprintf(stderr, "# Error: No parameter file specified!\n");
      exit(0);
    }
    
    strcpy(dnest_options_file_line, "\0"); /* empty the string */
    if(my_optind+1 < argc)
    {
      strcpy(dnest_options_file_line, argv[my_optind+1]); /* copy input OPTIONS1D file */
    }
    
    strcpy(dnest_options_file_con, "\0"); /* empty the string */
    if(my_optind+2 < argc)
    {
      strcpy(dnest_options_file_con, argv[my_optind+2]); /* copy input OPTIONSCON file */
    }
  }
  MPI_Bcast(&flag_postprc, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_end, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_restart, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_decomp, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_examine, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_load_prior, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_para_name, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_postsample, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  
  if(flag_examine == 1)
  {
    if(thistask == roottask)
    {
      test_covariance();
    }
    flag_end = 1;
  }

  if(flag_end != 1)
  {
    begin_run();

    end_run();
  }

  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
    int ht, mt;
    double st;
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
    printf("Ends successfully.\n");
    printf("===============MICA2==================\n");
  }
  return 0;
}
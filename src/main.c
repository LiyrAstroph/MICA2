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

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  double t0=0.0, t1=0.0, dt;
  int opt;
  extern int optind, opterr, optopt;
  extern char *optarg;
  extern int getopt(int argc, char *const *argv, const char *options);

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

    opterr = 0; /* reset getopt. */
    optind = 0; /* reset getopt. */
    
    flag_postprc = 0;
    flag_end = 0;
    flag_restart = 0;
    flag_decomp = 0;
    flag_examine = 0;
    flag_load_prior = 0;
    flag_para_name = 0;
    flag_postsample = 0;

    /* MAC getopt and GNU  getopt seem not compatible */
#if defined(__APPLE__) && defined(__MACH__)
    while( (opt = getopt(argc-1, argv+1, "pvrdel:ns")) != -1)
#else
    while( (opt = getopt(argc, argv, "pvrdel:ns")) != -1)
#endif
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
          strcpy(prior_file, optarg);
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
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;
      }
    }

    if(optind < argc ) // parameter file is specified 
      strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
    else
    {
      fprintf(stderr, "# Error: No parameter file specified!\n");
      exit(0);
    }
    
    strcpy(dnest_options_file_line, "\0"); /* empty the string */
    if(optind+1 < argc)
    {
      strcpy(dnest_options_file_line, argv[optind+1]); /* copy input OPTIONS1D file */
    }
    
    strcpy(dnest_options_file_con, "\0"); /* empty the string */
    if(optind+2 < argc)
    {
      strcpy(dnest_options_file_con, argv[optind+2]); /* copy input OPTIONSCON file */
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
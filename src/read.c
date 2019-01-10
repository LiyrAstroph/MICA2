/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/* read parset from param file */
int read_parset()
{
  if(thistask == roottask)
  {
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    nt = 0;
    strcpy(tag[nt], "FileDir");
    addr[nt] = &parset.file_dir;
    id[nt++] = STRING;

    strcpy(tag[nt], "DataFile");
    addr[nt] = &parset.data_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "FlagUniformVarParams");
    addr[nt] = &parset.flag_uniform_var_params;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagUniformTranFuns");
    addr[nt] = &parset.flag_uniform_tranfuns;
    id[nt++] = INT;

    strcpy(tag[nt], "LagLimitLow");
    addr[nt] = &parset.lag_limit_low;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "LagLimitUpp");
    addr[nt] = &parset.lag_limit_upper;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "NumGaussianLow");
    addr[nt] = &parset.num_gaussian_low;
    id[nt++] = INT;

    strcpy(tag[nt], "NumGaussianUpp");
    addr[nt] = &parset.num_gaussian_upper;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagConSysErr");
    addr[nt] = &parset.flag_con_sys_err;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagLineSysErr");
    addr[nt] = &parset.flag_line_sys_err;
    id[nt++] = INT;

    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    /* default parset */
    parset.num_gaussian_low = parset.num_gaussian_upper = 1;
    parset.flag_uniform_tranfuns = parset.flag_uniform_var_params = 0;
    parset.lag_limit_low = 0.0;
    parset.lag_limit_upper = -1.0;
    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='#')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
        {
          j = i;
          tag[i][0] = 0;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(id[j])
        {
          case DOUBLE:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy(addr[j], buf2);
            break;
          case INT:
            *((int *)addr[j]) = (int) atof(buf2);
            break;
        }
      }
      else
      {
        fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                      parset.param_file, buf1);
        exit(0);
      }
    }
    fclose(fparam);

    if(parset.num_gaussian_low <= 0)
    {
      printf("NumGaussianLow should be larger than 0.\n");
      exit(0);
    }

    if(parset.num_gaussian_upper <= 0)
    {
      printf("NumGaussianUpp should be larger than 0.\n");
      exit(0);
    }

    if(parset.num_gaussian_low > parset.num_gaussian_upper)
    {
      printf("NumGaussianLow should be smaller than or equal to NumGaussianUpp.\n");
      exit(0);
    }

    parset.num_gaussian_diff = parset.num_gaussian_upper - parset.num_gaussian_low + 1;
  }

  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return 0;
}


/* read date file */
int read_data()
{
  FILE *fp;
  char buf[256], str[256], str2[256], *pstr;
  int i, j, k, np;
  double tcad;

  /* read number of data sets. */
  if(thistask == roottask)
  { 
    sprintf(buf, "%s/%s", parset.file_dir, parset.data_file);
    fp = fopen(buf, "r");
    fgets(buf, 256, fp);
    sscanf(buf, "# %d\n", &nset);
  }

  MPI_Bcast(&nset, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  dataset = (DATASET *)malloc(nset * sizeof(DATASET));

  /* read number of continuum and number of line sets. */
  if(thistask == roottask)
  {
    for(i=0; i<nset; i++)
    {
      fgets(buf, 256, fp);
      sscanf(buf, "%s %s\n", str, str2);
      if(str[0]!='#')
      {
        printf("# Error.\n");
      }
      
      pstr = str2;
      sscanf(pstr, "%d", &(dataset[i].con.n));
      pstr = strchr(pstr, ':');

      dataset[i].nlset = 0;
      do 
      {
        pstr++;
        dataset[i].nlset++;
        pstr = strchr(pstr, ':');
      }while(pstr!=NULL);
    }
  }

  /* allocate memory. */
  for(i=0; i<nset; i++)
  {
    MPI_Bcast(&dataset[i].con.n, 1, MPI_INT, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&dataset[i].nlset, 1, MPI_INT, roottask, MPI_COMM_WORLD);

    dataset[i].con.t = malloc(dataset[i].con.n * sizeof(double));
    dataset[i].con.f = malloc(dataset[i].con.n * sizeof(double));
    dataset[i].con.fe = malloc(dataset[i].con.n * sizeof(double));

    dataset[i].line = malloc(dataset[i].nlset * sizeof(LC));
  }

  /* read number of points in each line light curves. */
  if(thistask == roottask)
  {
    rewind(fp);
    fgets(buf, 256, fp);

    for(i=0; i<nset; i++)
    {
      fgets(buf, 256, fp);
      sscanf(buf, "%s %s\n", str, str2);
      
      pstr = str2;
      pstr = strchr(pstr, ':');
      pstr++;
      sscanf(pstr, "%d", &(dataset[i].line[0].n));
      printf("%d %d ", dataset[i].con.n,  dataset[i].line[0].n);
      for(j=1; j<dataset[i].nlset; j++)
      {
        pstr = strchr(pstr, ':');
        pstr++;
        sscanf(pstr, "%d", &(dataset[i].line[j].n));
        printf("%d ", dataset[i].line[j].n);
      }
      printf("\n");
    }
  }

  /* allocate memory. */
  for(i=0; i<nset; i++)
  {
    for(j=0; j<dataset[i].nlset; j++)
    {
      MPI_Bcast(&(dataset[i].line[j].n), 1, MPI_INT, roottask, MPI_COMM_WORLD);

      dataset[i].line[j].t = malloc(dataset[i].line[j].n * sizeof(double));
      dataset[i].line[j].f = malloc(dataset[i].line[j].n * sizeof(double));
      dataset[i].line[j].fe = malloc(dataset[i].line[j].n * sizeof(double));
    }
  }

  /* now read light curve data. */
  if(thistask == roottask)
  {
    for(i=0; i<nset; i++)
    {
      /* continuum */
      for(j=0; j<dataset[i].con.n; j++)
      {
        fscanf(fp, "%lf %lf %lf\n", &(dataset[i].con.t[j]), &(dataset[i].con.f[j]), &(dataset[i].con.fe[j]));
      }
      fscanf(fp, "\n");

      /* line */
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<dataset[i].line[j].n; k++)
        {
          fscanf(fp, "%lf %lf %lf\n", &(dataset[i].line[j].t[k]), &(dataset[i].line[j].f[k]), &(dataset[i].line[j].fe[k]));
        }
        fscanf(fp, "\n");
      }
    }
  }

  /* broadcast data */
  for(i=0; i<nset; i++)
  {
    MPI_Bcast(dataset[i].con.t, dataset[i].con.n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(dataset[i].con.f, dataset[i].con.n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(dataset[i].con.fe, dataset[i].con.n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

    for(j=0; j<dataset[i].nlset; j++)
    {
      MPI_Bcast(dataset[i].line[j].t, dataset[i].line[j].n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
      MPI_Bcast(dataset[i].line[j].f, dataset[i].line[j].n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
      MPI_Bcast(dataset[i].line[j].fe, dataset[i].line[j].n, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    }
  }

  scale_con_line();
  cal_mean_error();

  
  /* combine all light curves */
  alldata = malloc(nset * sizeof(LC));
  for(i=0; i<nset; i++)
  {
    alldata[i].n = dataset[i].con.n;
    for(j=0; j<dataset[i].nlset; j++)
      alldata[i].n += dataset[i].line[j].n;

    alldata[i].t = malloc(alldata[i].n * sizeof(double));
    alldata[i].f = malloc(alldata[i].n * sizeof(double));
    alldata[i].fe = malloc(alldata[i].n * sizeof(double));

    memcpy(alldata[i].t, dataset[i].con.t, dataset[i].con.n*sizeof(double));
    memcpy(alldata[i].f, dataset[i].con.f, dataset[i].con.n*sizeof(double));
    memcpy(alldata[i].fe, dataset[i].con.fe, dataset[i].con.n*sizeof(double));
    np = dataset[i].con.n;
    for(j=0; j<dataset[i].nlset; j++)
    {
      memcpy(&alldata[i].t[np], dataset[i].line[j].t, dataset[i].line[j].n*sizeof(double));
      memcpy(&alldata[i].f[np], dataset[i].line[j].f, dataset[i].line[j].n*sizeof(double));
      memcpy(&alldata[i].fe[np], dataset[i].line[j].fe, dataset[i].line[j].n*sizeof(double));

      np += dataset[i].line[j].n;
    }
  }

  ncon_max = 0;
  nline_max = 0;
  nall_max = 0;
  nlset_max = 0;
  for(i=0; i<nset; i++)
  {
    if(ncon_max < dataset[i].con.n)
      ncon_max = dataset[i].con.n;

    if(nall_max < alldata[i].n)
      nall_max = alldata[i].n;

    if(nlset_max < dataset[i].nlset)
      nlset_max = dataset[i].nlset;

    for(j=0; j<dataset[i].nlset; j++)
    {
      if(nline_max < dataset[i].line[j].n)
        nline_max = dataset[i].line[j].n;
    }
  }

  tspan_max = 0.0;
  for(i=0; i<nset; i++)
  {
    if(tspan_max < dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0])
      tspan_max = dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0];
  }
  if(parset.lag_limit_upper < 0.0)
  {
    parset.lag_limit_upper = tspan_max/2.0;
  }

  tcadence_con_min = tspan_max;
  tcadence_line_min = tspan_max;
  for(i=0; i<nset; i++)
  {

    tcad = (dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0])/(dataset[i].con.n-1);
    if(tcadence_con_min > tcad)
      tcadence_con_min = tcad;

    for(j=0; j<dataset[i].nlset; j++)
    {
      tcad = (dataset[i].line[j].t[dataset[i].line[j].n-1] - dataset[i].line[j].t[0])/(dataset[i].line[j].n-1);
      if(tcadence_line_min  > tcad)
        tcadence_line_min = tcad;
    }
  }

  /* test */
  if(thistask == roottask)
  {
    int n;
    for(i=0; i<nset; i++)
    {
      printf("%f %f %f\n", dataset[i].con.t[0], dataset[i].con.f[0], dataset[i].con.fe[0]);
      printf("%f %f %f\n", dataset[i].con.t[dataset[i].con.n-1], dataset[i].con.f[dataset[i].con.n-1], dataset[i].con.fe[dataset[i].con.n-1]);

      n = dataset[i].line[0].n;
      printf("%f %f %f\n", dataset[i].line[0].t[0], dataset[i].line[0].f[0], dataset[i].line[0].fe[0]);
      printf("%f %f %f\n", dataset[i].line[0].t[n-1], dataset[i].line[0].f[n-1], dataset[i].line[0].fe[n-1]);
      printf("\n");
    }
  }

  return 0;
}

void cal_mean_error()
{
  int i, j, k;
  double mean;

  for(i=0; i<nset; i++)
  {
    /* continuum */
    mean = 0.0;
    for(j=0; j<dataset[i].con.n; j++)
    {
      mean += dataset[i].con.fe[j];
    }
    mean /= dataset[i].con.n;
    dataset[i].con.error_mean=mean;

    /* line */
    for(j=0; j<dataset[i].nlset; j++)
    {
      mean = 0.0;
      for(k=0; k<dataset[i].line[j].n; k++)
      {
        mean += dataset[i].line[j].fe[k];
      }
      mean /= dataset[i].line[j].n;
      dataset[i].line[j].error_mean=mean;
    }
  }

  return;
}

void scale_con_line()
{
  int i, j, k;
  double mean;

  for(i=0; i<nset; i++)
  {
    /* continuum */
    mean = 0.0;
    for(j=0; j<dataset[i].con.n; j++)
    {
      mean += dataset[i].con.f[j];
    }
    mean /= dataset[i].con.n;
    dataset[i].con.scale=mean;
    for(j=0; j<dataset[i].con.n; j++)
    {
      dataset[i].con.f[j] /= mean;
      dataset[i].con.fe[j] /= mean;
    }

    /* line */
    for(j=0; j<dataset[i].nlset; j++)
    {
      mean = 0.0;

      for(k=0; k<dataset[i].line[j].n; k++)
      {
        mean += dataset[i].line[j].f[k];
      }
      mean /= dataset[i].line[j].n;
      dataset[i].line[j].scale=mean;

      for(k=0; k<dataset[i].line[j].n; k++)
      {
        dataset[i].line[j].f[k] /= mean;
        dataset[i].line[j].fe[k] /= mean;
      }
    }
  }

  return;
}

/*!
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
{
  FILE *fp;
  char buf[MICA_MAX_STR_LENGTH], buf1[MICA_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, MICA_MAX_STR_LENGTH, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  sscanf(buf, "%d", &parset.num_particles);
  fclose(fp);
}

/*!
 * get file name of posterior sample. 
 */
void get_posterior_sample_file(char *fname, char *samplefile)
{
  FILE *fp;
  char buf[MICA_MAX_STR_LENGTH], buf1[MICA_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, MICA_MAX_STR_LENGTH, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.new_level_interval);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.save_interval);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.thread_steps);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_levels);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.lambda);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.beta);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_saves);

  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_file);
  
  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_info_file);
  
  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.levels_file);
  
  fgets(buf, MICA_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sampler_state_file);
  
  fgets(buf, MICA_MAX_STR_LENGTH, fp);
  sscanf(buf, "%s", samplefile);
  fclose(fp);
}

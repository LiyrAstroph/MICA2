/* MICA2
 * A code for time lag measurement in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 22, 2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "allvars.h"
#include "proto.h"

void set_param_file(char *fname)
{
  strcpy(parset.param_file, fname);
  return;
}

void set_parset(PARSET *ps)
{
  memcpy(&parset, ps, sizeof(PARSET));
  return;
}

void get_parset(PARSET *ps)
{
  memcpy(ps, &parset, sizeof(PARSET));
  return;
}

/* read parset from param file */
int read_parset()
{
  if(thistask == roottask)
  {
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam=NULL;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];

    pardict = malloc(MAXTAGS * sizeof(PARDICT));

    nt = 0;
    
    strcpy(pardict[nt].tag, "FileDir");
    pardict[nt].addr = &parset.file_dir;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "DataFile");
    pardict[nt].addr = &parset.data_file;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;
    
    strcpy(pardict[nt].tag, "TypeModel");
    pardict[nt].addr = &parset.model;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "MaxNumberSaves");
    pardict[nt].addr = &parset.max_num_saves;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagUniformVarParams");
    pardict[nt].addr = &parset.flag_uniform_var_params;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagUniformTranFuns");
    pardict[nt].addr = &parset.flag_uniform_tranfuns;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagLongtermTrend");
    pardict[nt].addr = &parset.flag_trend;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagNegativeResp");
    pardict[nt].addr = &parset.flag_negative_resp;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "LagLimitLow");
    pardict[nt].addr = &parset.lag_limit_low;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "LagLimitUpp");
    pardict[nt].addr = &parset.lag_limit_upper;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "WidthLimitLow");
    pardict[nt].addr = &parset.width_limit_low;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "WidthLimitUpp");
    pardict[nt].addr = &parset.width_limit_upper;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;
    
    /*===========================================*/
    /* to compatible with older version */
    strcpy(pardict[nt].tag, "NumGaussianLow");
    pardict[nt].addr = &parset.num_gaussian_low;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "NumGaussianUpp");
    pardict[nt].addr = &parset.num_gaussian_upper;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;
    /* to compatible with older version */
    /*===========================================*/

    strcpy(pardict[nt].tag, "NumCompLow");
    pardict[nt].addr = &parset.num_gaussian_low;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "NumCompUpp");
    pardict[nt].addr = &parset.num_gaussian_upper;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagConSysErr");
    pardict[nt].addr = &parset.flag_con_sys_err;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagLineSysErr");
    pardict[nt].addr = &parset.flag_line_sys_err;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "TypeLagPrior");
    pardict[nt].addr = &parset.type_lag_prior;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "StrLagPrior");
    pardict[nt].addr = &parset.str_lag_prior;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "StrWidthPrior");
    pardict[nt].addr = &parset.str_width_prior;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "StrRatioPrior");
    pardict[nt].addr = &parset.str_ratio_prior;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "TypeTF");
    pardict[nt].addr = &parset.type_tf;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "StrTypeTFMix");
    pardict[nt].addr = &parset.str_type_tf_mix;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "FlagLagPositivity");
    pardict[nt].addr = &parset.flag_lag_posivity;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "FlagGap");
    pardict[nt].addr = &parset.flag_gap;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "StrGapPrior");
    pardict[nt].addr = &parset.str_gap_prior;
    pardict[nt].isset = 0;
    pardict[nt++].id = STRING;

    strcpy(pardict[nt].tag, "NumberParticles");
    pardict[nt].addr = &parset.num_particles;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "NewLevelIntervalFactor");
    pardict[nt].addr = &parset.new_level_interval_factor;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;
  
    strcpy(pardict[nt].tag, "SaveIntervalFactor");
    pardict[nt].addr = &parset.save_interval_factor;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;
  
    strcpy(pardict[nt].tag, "ThreadStepsFactor");
    pardict[nt].addr = &parset.thread_steps_factor;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;
  
    strcpy(pardict[nt].tag, "MaxNumberLevels");
    pardict[nt].addr = &parset.max_num_levels;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;
  
    strcpy(pardict[nt].tag, "BacktrackingLength");
    pardict[nt].addr = &parset.lam;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;
   
    strcpy(pardict[nt].tag, "StrengthEqualPush");
    pardict[nt].addr = &parset.beta;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "PTol");
    pardict[nt].addr = &parset.max_ptol;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "NumPointRec");
    pardict[nt].addr = &parset.nd_rec;
    pardict[nt].isset = 0;
    pardict[nt++].id = INT;

    strcpy(pardict[nt].tag, "TimeRecLowExt");
    pardict[nt].addr = &parset.trec_low_ext;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    strcpy(pardict[nt].tag, "TimeRecUppExt");
    pardict[nt].addr = &parset.trec_upp_ext;
    pardict[nt].isset = 0;
    pardict[nt++].id = DOUBLE;

    num_pardict = nt;
    
    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    /* default parset */
    strcpy(parset.file_dir,"\0"); /* empty file dir */
    parset.num_gaussian_low = parset.num_gaussian_upper = 1;
    parset.flag_uniform_tranfuns = parset.flag_uniform_var_params = 0;
    parset.lag_limit_low = 0.0;
    parset.lag_limit_upper = -1.0;
    parset.type_lag_prior = 0;
    parset.width_limit_low = DBL_MAX;
    parset.width_limit_upper = DBL_MAX;
    parset.flag_trend = 0;
    parset.flag_negative_resp = 0;
    parset.type_tf = 0;
    parset.flag_lag_posivity = 0;
    parset.num_gaussian_low = 1;
    parset.num_gaussian_upper = 1;
    parset.flag_gap = 0;
    parset.nd_rec = 200;
    parset.trec_low_ext = 0.0;
    parset.trec_upp_ext = 0.0;
    strcpy(parset.str_lag_prior,"");
    strcpy(parset.str_ratio_prior,"");
    strcpy(parset.str_width_prior,"");
    strcpy(parset.str_gap_prior,"");
    strcpy(parset.str_type_tf_mix,"");
    /*cdnest options */
    parset.num_particles = 1;
    parset.max_num_saves = 2000;
    parset.new_level_interval_factor = 1;
    parset.save_interval_factor = parset.new_level_interval_factor;
    parset.thread_steps_factor = 1;
    parset.lam = 10.0;
    parset.beta = 100.0;
    parset.max_num_levels = 0;
    parset.max_ptol = 0.1;

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='#')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, pardict[i].tag) == 0 && pardict[i].isset == 0)
        {
          j = i;
          pardict[i].isset = 1;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(pardict[j].id)
        {
          case DOUBLE:
            *((double *) pardict[j].addr) = atof(buf2);
            break;
          case STRING:
            strcpy(pardict[j].addr, buf2);
            break;
          case INT:
            *((int *)pardict[j].addr) = (int) atof(buf2);
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

    if((parset.model<gmodel) || (parset.model>nmap))
    {
      printf("TypeModel should be 0-4.\n 0: general model; 1: pmap; 2: vmap; 3: mmap; 4: nmap.\n");
      exit(0);
    }

    /* do not check options for nmap, because they are not used */
    if(parset.model != nmap)
    {
      if(parset.model == 0)
      {
        if(parset.lag_limit_upper <= parset.lag_limit_low)
        {
          printf("LagLimitUpp should be larger than LagLimitLow!\n");
          exit(0);
        }
      }

      if(parset.width_limit_low < DBL_MAX)
      {
        parset.width_limit_low_isset = 1;
      }
      if(parset.width_limit_upper < DBL_MAX)
      {
        parset.width_limit_upper_isset = 1;
      }

      if(parset.width_limit_low_isset == 1) /* is set */
      {
        if(parset.width_limit_low <= 0.0)
        {
          printf("WidthLimitLow soubld be positive!\n");
          exit(0);
        }
      }
      if(parset.width_limit_upper_isset == 1) /* is set */
      {
        if(parset.width_limit_upper <= 0.0)
        {
          printf("WidthLimitUpp soubld be positive!\n");
          exit(0);
        }
      }
      if(parset.width_limit_upper_isset == 1 && parset.width_limit_low_isset == 1)
      {
        if(parset.width_limit_upper <= parset.width_limit_low)
        {
          printf("WidthLimitUpp should be larger than WidthLimitLow!\n");
          exit(0);
        }
      }

      /* pmap mode, num_gaussian >= 2 */
      if(parset.model == pmap)
      {
        if(parset.num_gaussian_low < 2)
          parset.num_gaussian_low = 2;

        if(parset.num_gaussian_upper < 2)
          parset.num_gaussian_upper = 2;
        
        parset.type_lag_prior = 4;
        if (parset.str_lag_prior[0] == '\0')
        {
          printf("StrLagPrior should be set for pmap mode!\n");
          exit(0);
        }
      }

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

      if(parset.num_gaussian_upper == 1)
      {
        parset.type_lag_prior = 0;
      }
      
      /* when loading prior, only support a single value of num_gaussian */
      if(flag_load_prior == 1 && parset.num_gaussian_diff > 1)
      {
        printf("For loading prior (-l), only support NumGausssianUpp=NumGaussianLow.\n");
        exit(0);
      }

      if(parset.nd_rec < 0)
      {
        printf("NumPointRec should be larger than 0.\n");
        exit(0);
      }

      if(parset.type_lag_prior < 0 || parset.type_lag_prior > 4)
      {
        printf("Incorrect TypeLagPrior, should be 0-4.\n");
        exit(0);
      }

      if(parset.type_lag_prior >= 2 && parset.type_lag_prior <=3 && parset.num_gaussian_low < 3)
      {
        printf("For prior type 2 or 3, better to use more Gaussians (>=3)!\n");
        exit(0);
      }

      if(parset.flag_trend < 0)
      {
        printf("Incorrect FlagLongtermTrend, should be equal to or larger than 0.\n");
        exit(0);
      }
      if(parset.flag_trend > 2)
      {
        printf("FlagLongtermTrend seems too large.\n");
        exit(0);
      }
      if(parset.flag_lag_posivity != 0 && parset.lag_limit_low < 0.0)
      {
        printf("LagLimitLow shoule be >=0 when enabling FlagLagPositivity.\n");
        exit(0);
      }
      if(parset.flag_negative_resp < 0 && parset.flag_negative_resp > 1)
      {
        printf("FlagNegativeResp shoule be 0 or 1.\n");
        exit(0);
      }
      if(parset.flag_negative_resp == 1 && parset.model != 0)
      {
        printf("FlagNegativeResp only works when TypeModel = 0.\n");
        exit(0);
      }
      if(parset.flag_negative_resp == 1 && parset.type_lag_prior == 0)
      {
        /* make sure the lag prior type > 0 */
        parset.type_lag_prior = 1;
      }
      if(parset.flag_lag_posivity == 1 && (parset.type_tf > 1) && (parset.type_tf <= 3))
      {
        printf("For gamma and exp tf, no need to use FlagLagPositivity, instead, set LagLimitLow > 0.\n");
        exit(0);
      }
      if(parset.flag_uniform_tranfuns != 0 && parset.model != 0)
      {
        printf("For uniform transfuns, typemodel should be 0 (general model).\n");
        exit(0);
      }
      if(parset.model == mmap)
      {
        /* presently, do not check posivity for mmap, may be adjusted late on */
        parset.flag_lag_posivity = 0; 
        if(parset.str_type_tf_mix[0] == '\0')
        {
          printf("For mmap mode, StrTypeTFMix should be set!\n");
          exit(-1);
        }

        int len = strlen(parset.str_type_tf_mix);
        if(len < 2)
        {
          printf("StrTypeTFMix (='%s') should contain at least two numbers!\n", parset.str_type_tf_mix);
          exit(-1);
        }

        /* check whether StrTypeTFMix is valid */
        for(i=0; i<len; i++)
        {
          if(parset.str_type_tf_mix[i] < '0' || parset.str_type_tf_mix[i] > '3')
          {
            printf("StrTypeTFMix (='%s') should only contain numbers 0-3!\n", parset.str_type_tf_mix);
            exit(-1);
          }
        }
        
        /* fix the number of components */
        parset.num_gaussian_low = len;
        parset.num_gaussian_upper = len;
      }
    }
  }

  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return 0;
}


/* read date file */
int read_data()
{
  FILE *fp=NULL;
  char buf[256], *pstr;
  int i, j, k, np;
  double tcad, tspan;

  /* read number of data sets. */
  if(parset.model != nmap)
  {
    if(thistask == roottask)
    { 
      sprintf(buf, "%s/%s", parset.file_dir, parset.data_file);
      fp = fopen(buf, "r");
      if(fp == NULL)
      {
        printf("Cannot open file %s.\n", buf);
        exit(0);
      }
      fgets(buf, 256, fp);
      if(buf[0]!='#')
      {
        printf("# Error in input data.\n");
      }
      sscanf(buf+1, "%d\n", &nset);
    }

    MPI_Bcast(&nset, 1, MPI_INT, roottask, MPI_COMM_WORLD);

    dataset = (DATASET *)malloc(nset * sizeof(DATASET));

    /* read number of continuum and number of line sets. */
    if(thistask == roottask)
    {
      for(i=0; i<nset; i++)
      {
        fgets(buf, 256, fp);
        if(buf[0]!='#')
        {
          printf("# Error in input data.\n");
        }
        
        pstr = buf+1;
        sscanf(pstr, "%d", &(dataset[i].con.n));
        pstr = strchr(pstr, ':');
        if(pstr == NULL) 
        {
          dataset[i].nlset = 0;
          continue;
        }

        dataset[i].nlset = 0;
        do 
        {
          pstr++;
          
          /* remove blanks */
          while(*pstr == ' ' || *pstr == '\t') pstr++;
          if(*pstr == '\n' || *pstr == '\r' || *pstr == '\0')
            *pstr = '\0';

          if(*pstr=='\0')
            np = 0;
          else
            sscanf(pstr, "%d", &np);
            
          if(np > 0) /* only account those lines with points*/
          {
            dataset[i].nlset++;
          }
          pstr = strchr(pstr, ':');
        }while(pstr!=NULL);

        if(dataset[i].con.n >= 1000)
        {
          printf("The number of points in %d-th dataset is %d, a bit large!\n", i, dataset[i].con.n);
          printf("Better to rebin the light curve to reduce the number of points.\n"
                "This will improve the computational speed [~O(N^3)]!\n");
          // exit(-1);
        }
        
        /* for vmap, zero continuum point */
        if(parset.model == vmap)
        {
          if(dataset[i].con.n > 0)
          {
            printf("vmap mode, the number of continuum points in %d-th dataset should be zero!\n", i);
            exit(-1);
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* allocate memory. */
    for(i=0; i<nset; i++)
    {
      MPI_Bcast(&dataset[i].con.n, 1, MPI_INT, roottask, MPI_COMM_WORLD);
      MPI_Bcast(&dataset[i].nlset, 1, MPI_INT, roottask, MPI_COMM_WORLD);

      dataset[i].con.t  = (double *)malloc(dataset[i].con.n * sizeof(double));
      dataset[i].con.f  = (double *)malloc(dataset[i].con.n * sizeof(double));
      dataset[i].con.fe = (double *)malloc(dataset[i].con.n * sizeof(double));

      dataset[i].line   = (LC *)malloc(dataset[i].nlset * sizeof(LC));
    }

    /* read number of points in each line light curves. */
    if(thistask == roottask)
    {
      rewind(fp);
      fgets(buf, 256, fp);

      for(i=0; i<nset; i++)
      {
        printf("set %d, # of points\n", i);
        fgets(buf, 256, fp);
        pstr = buf+1; /* skip the first # */

        printf("%d ", dataset[i].con.n);
        for(j=0; j<dataset[i].nlset; j++)
        {
          pstr = strchr(pstr, ':');
          pstr++;
          sscanf(pstr, "%d", &np);
          if(np > 0)
          {
            dataset[i].line[j].n = np;
            printf("%d ", dataset[i].line[j].n);
          }
        }
        printf("\n");
      }
    }
  }
  else /* for nmap, treat all light curves as continuum */
  {
    int nset_temp, iset;
    if(thistask == roottask)
    { 
      sprintf(buf, "%s/%s", parset.file_dir, parset.data_file);
      fp = fopen(buf, "r");
      if(fp == NULL)
      {
        printf("Cannot open file %s.\n", buf);
        exit(0);
      }
      fgets(buf, 256, fp); 
      sscanf(buf+1, "%d\n", &nset_temp);
      
      /* first read the number of light curves */
      nset = 0;
      for(i=0; i<nset_temp; i++)
      {
        fgets(buf, 256, fp);
        if(buf[0]!='#')
        {
          printf("# Error in input data.\n");
        }

        pstr = buf+1; /* skip the first # */
        nset++;
        pstr = strchr(pstr, ':');
        if(pstr == NULL) 
        {
          dataset[i].nlset = 0;
          continue;
        }

        do 
        {
          pstr++;

          /* remove blanks */
          while(*pstr == ' ' || *pstr == '\t') pstr++;
          if(*pstr == '\n' || *pstr == '\r' || *pstr == '\0')
            *pstr = '\0';
          
          if(*pstr=='\0')  
            np = 0;
          else
            sscanf(pstr, "%d", &np);

          if(np > 0) /* only account those lines with points */
          {
            nset++;
          }
          pstr = strchr(pstr, ':');
        }while(pstr!=NULL);
      }
    }
    
    MPI_Bcast(&nset, 1, MPI_INT, roottask, MPI_COMM_WORLD);
    dataset = (DATASET *)malloc(nset * sizeof(DATASET));
    
    if(thistask == roottask)
    {
      /* now read the number of points in each light curve */
      rewind(fp);
      fgets(buf, 256, fp); /* the first line is no longer used for nmap */
      
      iset = 0;
      for(i=0; i<nset_temp; i++)
      {
        fgets(buf, 256, fp);
        if(buf[0]!='#')
        {
          printf("# Error in input data.\n");
        }

        pstr = buf+1; /* skip the first # */
        sscanf(pstr, "%d", &(dataset[iset].con.n));
        printf("Set %d, %d points\n", iset, dataset[iset].con.n);
        iset++;
        pstr = strchr(pstr, ':');
        if(pstr == NULL) 
        {
          continue;
        }

        do
        {
          pstr++;

          /* remove blanks */
          while(*pstr == ' ' || *pstr == '\t') pstr++;
          if(*pstr == '\n' || *pstr == '\r' || *pstr == '\0')
            *pstr = '\0';
          
          if(*pstr=='\0')  
            np = 0;
          else
            sscanf(pstr, "%d", &np);
            
          if(np > 0) /* only account those lines with points */
          {
            dataset[iset].con.n = np;
            printf("Set %d, %d points\n", iset, dataset[iset].con.n);
            iset++;
          }
          pstr = strchr(pstr, ':');
        }while(pstr!=NULL);
      }
      
      for(i=0; i<nset; i++)
      {
        dataset[i].nlset = 0; /* no line set */
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i<nset; i++)
    {
      MPI_Bcast(&dataset[i].con.n, 1, MPI_INT, roottask, MPI_COMM_WORLD);
      MPI_Bcast(&dataset[i].nlset, 1, MPI_INT, roottask, MPI_COMM_WORLD);

      dataset[i].con.t  = (double *)malloc(dataset[i].con.n * sizeof(double));
      dataset[i].con.f  = (double *)malloc(dataset[i].con.n * sizeof(double));
      dataset[i].con.fe = (double *)malloc(dataset[i].con.n * sizeof(double));

      dataset[i].line   = (LC *)malloc(dataset[i].nlset * sizeof(LC));
      for(j=0; j<dataset[i].nlset; j++)
      {
        dataset[i].line[j].n = 0; /* no line set */
      }
    }
  }
  
  /* allocate memory. */
  for(i=0; i<nset; i++)
  {    
    for(j=0; j<dataset[i].nlset; j++)
    {
      MPI_Bcast(&(dataset[i].line[j].n), 1, MPI_INT, roottask, MPI_COMM_WORLD);

      dataset[i].line[j].t  = (double *)malloc(dataset[i].line[j].n * sizeof(double));
      dataset[i].line[j].f  = (double *)malloc(dataset[i].line[j].n * sizeof(double));
      dataset[i].line[j].fe = (double *)malloc(dataset[i].line[j].n * sizeof(double));
      if(dataset[i].line[j].n == 0)
      {
        parset.model = nmap;
      }
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
        if(feof(fp) != 0)
        {
          printf("error in reading the data file %s.\n", parset.data_file);
          exit(0);
        }
        np = fscanf(fp, "%lf %lf %lf\n", &(dataset[i].con.t[j]), &(dataset[i].con.f[j]), &(dataset[i].con.fe[j]));
        if(np < 3)
        {
          printf("Error in reading %d-th row of continuum of %d-th dataset.\n", j, i);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      /* check whether time is sorted */
      if(check_time_sorted(dataset[i].con.t, dataset[i].con.n)==1)
      {
        printf("--- continuum of %d-th set.\n", i);
        exit(0);
      }

      /* line */
      for(j=0; j<dataset[i].nlset; j++)
      {
        for(k=0; k<dataset[i].line[j].n; k++)
        {
          if(feof(fp) != 0)
          {
            printf("error in reading the data file %s.\n", parset.data_file);
            exit(0);
          }  
          np = fscanf(fp, "%lf %lf %lf\n", &(dataset[i].line[j].t[k]), &(dataset[i].line[j].f[k]), &(dataset[i].line[j].fe[k]));
          if(np < 3)
          {
            printf("Error in reading %d-th row of %d-th line of %d-th dataset.\n", j, k, i);
            exit(0);
          }
        }
        fscanf(fp, "\n");

        /* check whether time is sorted */
        if(check_time_sorted(dataset[i].line[j].t, dataset[i].line[j].n) == 1)
        {
          printf("--- %d-th line of %d-th set.\n", j, i);
          exit(0);
        }
      }
    }
    fclose(fp);
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

  /* save max and min flux of each dataset */
  flux_minmax = (double **)malloc(nset * sizeof(double *));
  for(i=0; i<nset; i++)
  {
    flux_minmax[i] = (double *)malloc((1+dataset[i].nlset)*2 * sizeof(double));
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
  nrec_max = 0;
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

    nrec_max = fmax(nrec_max, fmax(ncon_max, nline_max) * (1 + dataset[i].nlset));
  }

  tspan_max = 0.0;
  for(i=0; i<nset; i++)
  {
    /* note that continuum might be empty */
    if(dataset[i].con.n > 0)
    {
      tspan = dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0];
      if(tspan_max < tspan)
        tspan_max = tspan;
    }
    
    for(j=0; j<dataset[i].nlset; j++)
    {
      if(dataset[i].line[j].n == 0)
        continue;

      tspan = dataset[i].line[j].t[dataset[i].line[j].n-1] - dataset[i].line[j].t[0];
      if(tspan_max < tspan)
        tspan_max = tspan;
    }
  }
  if(parset.lag_limit_upper < 0.0)
  {
    parset.lag_limit_upper = tspan_max/2.0;
  }

  tcadence_con_min = tspan_max;
  tcadence_line_min = tspan_max;
  for(i=0; i<nset; i++)
  {
    /* note that continuum might be empty */
    if(dataset[i].con.n > 0)
    {
      tcad = (dataset[i].con.t[dataset[i].con.n-1] - dataset[i].con.t[0])/(dataset[i].con.n-1);
      if(tcadence_con_min > tcad)
        tcadence_con_min = tcad;
    }

    for(j=0; j<dataset[i].nlset; j++)
    {
      if(dataset[i].line[j].n == 0)
        continue;

      tcad = (dataset[i].line[j].t[dataset[i].line[j].n-1] - dataset[i].line[j].t[0])/(dataset[i].line[j].n-1);
      if(tcadence_line_min  > tcad)
        tcadence_line_min = tcad;
    }
  }
  tcadence_min = fmin(tcadence_con_min, tcadence_line_min);

  /* test */
  if(thistask == roottask)
  {
    int n;
    printf("#=======================================================\n");
    for(i=0; i<nset; i++)
    {
      printf("set %d, starting and ending points\n", i);
      printf("%-4d %f %f %f\n", 0, dataset[i].con.t[0], dataset[i].con.f[0], dataset[i].con.fe[0]);
      printf("%-4d %f %f %f\n", dataset[i].con.n-1, dataset[i].con.t[dataset[i].con.n-1], dataset[i].con.f[dataset[i].con.n-1], dataset[i].con.fe[dataset[i].con.n-1]);
      for(j=0; j<dataset[i].nlset; j++)
      {
        n = dataset[i].line[j].n;
        if(n == 0)
          continue;

        printf("%-4d %f %f %f\n", 0, dataset[i].line[j].t[0], dataset[i].line[j].f[0], dataset[i].line[j].fe[0]);
        printf("%-4d %f %f %f\n", dataset[i].line[j].n-1, dataset[i].line[j].t[n-1], dataset[i].line[j].f[n-1], dataset[i].line[j].fe[n-1]);
      }
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
    if(parset.model != vmap)
    {
      mean = 0.0;
      for(j=0; j<dataset[i].con.n; j++)
      {
        mean += dataset[i].con.fe[j];
      }
      mean /= dataset[i].con.n;
      dataset[i].con.error_mean=mean;
    }
    else 
    {
      dataset[i].con.error_mean=0.0;
    }

    /* line */
    for(j=0; j<dataset[i].nlset; j++)
    {
      /* no data point */
      if(dataset[i].line[j].n == 0)
      {
        dataset[i].line[j].error_mean=0.0;
      }
      else 
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
  }

  return;
}

void scale_con_line()
{
  int i, j, k;
  double mean, Rmax, fdmax, fdmin, scale;  /* mean flux, difference between max and min fluxes*/
  int flag_norm = 0; /* by default use mean to do normalization  */

  /* first determine whether using mean or Rmax to do normalizations */
  if(thistask == roottask)
  {
    for(i=0; i<nset; i++)
    {
      if(parset.model != vmap )
      {
        /* continuum */
        mean = 0.0;
        fdmax = fdmin = dataset[i].con.f[0]; 
        for(j=0; j<dataset[i].con.n; j++)
        {
          mean += dataset[i].con.f[j];

          if(fdmax<dataset[i].con.f[j]) fdmax = dataset[i].con.f[j];
          if(fdmin>dataset[i].con.f[j]) fdmin = dataset[i].con.f[j];
        }
        mean /= dataset[i].con.n;
        Rmax = (fdmax-fdmin)/2.0;

        /* check if the mean is positive*/
        if(Rmax == 0.0)
        {
          printf("Error: the continuum of %d-th set is constant, without variability. \n", i);
          exit(0);
        }
        else if(mean < Rmax/100.0)
        {
          flag_norm = 1;
        }
        else
        {
          flag_norm = 0;
        }
      }

      /* line */
      for(j=0; j<dataset[i].nlset; j++)
      {
        /* no data point, continuue */
        if(dataset[i].line[j].n == 0)
          continue;

        mean = 0.0;
        fdmax = fdmin = dataset[i].line[j].f[0];
        for(k=0; k<dataset[i].line[j].n; k++)
        {
          mean += dataset[i].line[j].f[k];

          if(fdmax<dataset[i].line[j].f[k]) fdmax = dataset[i].line[j].f[k];
          if(fdmin>dataset[i].line[j].f[k]) fdmin = dataset[i].line[j].f[k];
        }
        mean /= dataset[i].line[j].n;
        Rmax = (fdmax-fdmin)/2.0;

        /* check if the mean is positive*/
        if(Rmax == 0.0)
        {
          printf("Error: the %d-th line of %d-th set is constant, without variability. \n", j, i);
          exit(0);
        }
        else if(mean < Rmax/100.0 ) /* use Rmax to do normalization */
        {
          flag_norm = 1;
        }
        else if (flag_norm != 1) /* use mean to do normalization, not yet set to 1 */
        {
          flag_norm = 0;
        }
      }
    }
  }
  
  /* now broadcast the value */
  MPI_Bcast(&flag_norm, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  /* do the normalization */
  for(i=0; i<nset; i++)
  {
    if(parset.model != vmap )
    {
      /* continuum */
      mean = 0.0;
      fdmax = fdmin = dataset[i].con.f[0]; 
      for(j=0; j<dataset[i].con.n; j++)
      {
        mean += dataset[i].con.f[j];

        if(fdmax<dataset[i].con.f[j]) fdmax = dataset[i].con.f[j];
        if(fdmin>dataset[i].con.f[j]) fdmin = dataset[i].con.f[j];
      }
      mean /= dataset[i].con.n;
      Rmax = (fdmax-fdmin)/2.0;

      if(flag_norm == 1) /* use Rmax to do normalization */
      {
        scale = Rmax;
        if(thistask == 0)
        {
          printf("Calculate Rmax=fmax-fmin (=%f) of continuum of %d-th set!\n", Rmax, i);
        }
      }
      else /* use mean to do normalization */
      {
        scale = mean;
        if(thistask == 0)
        {
          printf("Calculate the mean (=%f) of continuum of %d-th set!\n", mean, i);
        }
      }

      if(parset.flag_uniform_tranfuns == 1 || parset.flag_uniform_var_params == 1)
      {
        scale = 1.0;  /* no scale */
      }

      if(thistask == 0)
      {
        printf("Use %f to normalize the continuum light curve of %d-th set!\n", scale, i);
      }

      dataset[i].con.scale=scale;
      for(j=0; j<dataset[i].con.n; j++)
      {
        dataset[i].con.f[j] /= scale;
        dataset[i].con.fe[j] /= scale;
      }

      /* save flux min and max */
      flux_minmax[i][0*2 + 0] = (fdmin+fdmax)/2/scale;
      flux_minmax[i][0*2 + 1] = (fdmax-fdmin)/2/scale;
    }
    else /* for vmap, no continuum data points */
    {
      dataset[i].con.scale = 1.0;
    }

    /* line */
    for(j=0; j<dataset[i].nlset; j++)
    {
      if(dataset[i].line[j].n == 0)
      {
        dataset[i].line[j].scale=1.0;
        continue;
      }

      mean = 0.0;
      fdmax = fdmin = dataset[i].line[j].f[0];
      for(k=0; k<dataset[i].line[j].n; k++)
      {
        mean += dataset[i].line[j].f[k];

        if(fdmax<dataset[i].line[j].f[k]) fdmax = dataset[i].line[j].f[k];
        if(fdmin>dataset[i].line[j].f[k]) fdmin = dataset[i].line[j].f[k];
      }
      mean /= dataset[i].line[j].n;
      Rmax = (fdmax-fdmin)/2.0;

      if(flag_norm == 1) /* use Rmax to do normalization */
      {
        scale = Rmax;
        if(thistask == 0)
        {
          printf("Calculate Rmax=fmax-fmin (=%f) of the %d-th line of %d-th set!\n", Rmax, j, i);
        }
      }
      else /* use mean to do normalization */
      {
        scale = mean;
        if(thistask == 0)
        {
          printf("Calculate the mean (=%f) of the %d-th line of %d-th set !\n", mean, j, i);
        }
      }

      if(parset.flag_uniform_tranfuns == 1 || parset.flag_uniform_var_params == 1)
      {
        scale = 1.0;  /* no scale */
      }

      if(thistask == 0)
      {
        printf("Use %f to normalize the %d-th line of %d-th set !\n", scale, j, i);
      }

      dataset[i].line[j].scale=scale;

      for(k=0; k<dataset[i].line[j].n; k++)
      {
        dataset[i].line[j].f[k] /= scale;
        dataset[i].line[j].fe[k] /= scale;
      }

      /* save flux min and max */
      flux_minmax[i][(1+j)*2 + 0] = (fdmin+fdmax)/2/scale;
      flux_minmax[i][(1+j)*2 + 1] = (fdmax-fdmin)/2/scale;
    }
  }

  return;
}

/*!
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
{
  FILE *fp=NULL;
  char buf[MICA_MAX_STR_LENGTH], buf1[MICA_MAX_STR_LENGTH], buf2[MICA_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  /* default number particles */
  parset.num_particles = 1;

  while(!feof(fp))
  {
    fgets(buf, MICA_MAX_STR_LENGTH, fp);
    if(buf[0] == '#')
      continue;
    if(sscanf(buf, "%s%s", buf1, buf2) < 1)  /* a blank line */
      continue;
    if(sscanf(buf, "%s%s", buf1, buf2) < 2)
    {
      fprintf(stderr, "Error in geting number of particles.\n"
                      "Usually due to incorrect options.\n");
      exit(0);
    }
    if(strcmp(buf1, "NumberParticles") == 0)
    {
      parset.num_particles = atoi(buf2);
      break;
    }
  }

  fclose(fp);
}

int cmpfunc(const void * a, const void * b) 
{
  if( *(double*)a > *(double*)b ) return 1;
  if( *(double*)a < *(double*)b ) return -1;
  return 0;
}

/*!
 * get seasonal gap from continuum light curve data
 *
 * note that tcon must be monotonically increasing
 */
double get_seasonal_gap(const double *tcon, int ncon)
{
  int i, ny;
  double *dt;
  double span, gap;
  
  /* number of years */
  span = tcon[ncon-1] - tcon[0];
  ny = (int)ceil(span/YEAR_DAY);
  //printf("%d\n", ny);
  
  gap = 0.0;
  if(ny > 1)
  {
    dt = malloc((ncon-1)*sizeof(double));

    /* sampling interval */
    for(i=0; i<ncon-1; i++)
    {
      dt[i] = tcon[i+1] - tcon[i];
    }

    qsort(dt, ncon-1, sizeof(double), cmpfunc);

    /* for ny years, so there must be (ny-1) gaps */
    for(i=ncon-2; i>ncon-2-(ny-1); i--)
    {
      gap += dt[i];
    }
    gap /= (ny-1);

    free(dt);
  }
  else
  {
    gap = tcadence_con_min; 
  }

  return gap;
}

void get_seasonal_gap_allset()
{
  int i;

  if(strlen(parset.str_gap_prior) == 0)  /* default */
  {
    for(i=0; i<nset; i++)
    {
      gap_center[i] = YEAR_DAY/2.0;
      gap_width[i] = get_seasonal_gap(dataset[i].con.t, dataset[i].con.n);
      if(thistask == roottask)
        printf("Gap of set %d: %f %f\n", i, gap_center[i], gap_width[i]);
      gap_width[i] /= 2.0;
    }
  }
  else /* input from str_gap_prior */
  {
    char *pstr = parset.str_gap_prior;
    int j;

    pstr += 1;
    j = 0;
    for(i=0; i<nset-1; i++)
    {
      if(sscanf(pstr, "%lf:%lf", &gap_center[j], &gap_width[j]) < 2)
      {
        if(thistask == 0)
          printf("No enough gap priors.\n");
        exit(0);
      }
      if(thistask == roottask)
        printf("Gap of set %d: %f %f\n", i, gap_center[i], gap_width[i]);
      gap_width[i] /= 2.0;

      j++;

      pstr = strchr(pstr, ':'); /* values are separated by ":" */
      if(pstr!=NULL)
      {
        pstr++;
      }
      else
      {
        if(thistask == 0)
          printf("No enough gap priors.\n");
        exit(0);
      }
    }

    if(sscanf(pstr, "%lf:%lf", &gap_center[j], &gap_width[j]) < 2)
    {
      if(thistask == 0)
        printf("No enough gap priors.\n");
      exit(0);
    }
    if(thistask == roottask)
      printf("Gap of set %d: %f %f\n", j, gap_center[j], gap_width[j]);
    gap_width[j] /= 2.0;
  }
  return;
}

/*!
 * check whether data is sorted in time.
 * 
 */
int check_time_sorted(double *time_series, int n)
{
  int i;
  double dt;

  for(i=1; i<n; i++)
  {
    dt = time_series[i] - time_series[i-1];
    if(dt < 0.0)
    {
      printf("light curve data are not in ascending order in time.\n");
      return 1;
    }
  }
  return 0;
}
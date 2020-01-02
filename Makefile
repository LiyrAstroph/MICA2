


SHELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O2 -Wall -finline-functions
#OPTIMIZE += -DDebug

#------------target system---------
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"
#SYSTEM="TianheII"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
DNEST_INCL  = -I /home/liyropt/Projects/GIT/CDNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/CDNest -ldnest

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/sharefs/mbh/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/sharefs/mbh/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/sharefs/mbh/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/sharefs/mbh/user/liyanrong/soft/mpich3/include
LAPACK_INCL = -I/sharefs/mbh/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/sharefs/mbh/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran

DNEST_INCL  = -I /sharefs/mbh/user/liyanrong/GIT/DNest/
DNEST_LIBS  = -L /sharefs/mbh/user/liyanrong/GIT/DNest -ldnest
endif

ifeq ($(SYSTEM), "TianheII")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas -lm
MPICHLIB = -lmpich
MPIINCL  =
LAPACK_INCL = -I/HOME/ihep_yrli_1/BIGDATA/soft/lapack/include
LAPACK_LIBS = -L/HOME/ihep_yrli_1/BIGDATA/soft/lapack/lib -llapacke -llapack -lblas -lgfortran

DNEST_INCL  = -I /HOME/ihep_yrli_1/BIGDATA/soft/DNest/
DNEST_LIBS  = -L /HOME/ihep_yrli_1/BIGDATA/soft/DNest -ldnest
endif


EXEC     = mica2
SRC      = ./src
OBJS     = $(SRC)/main.o $(SRC)/allvars.o $(SRC)/system.o $(SRC)/run.o        \
           $(SRC)/dnest_con.o  $(SRC)/dnest_line.o $(SRC)/read.o              \
           $(SRC)/mc_con.o $(SRC)/init.o $(SRC)/mathfun.o                     \
           $(SRC)/mc_line.o  $(SRC)/error.o
  
INCL     = Makefile $(SRC)/allvars.h $(SRC)/proto.h  $(SRC)/dnest_con.h \
           $(SRC)/dnest_line.h 

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPICHINCL) $(DNEST_INCL) $(FFTW_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(MPICHLIB) $(DNEST_LIBS) $(FFTW_LIBS)

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)

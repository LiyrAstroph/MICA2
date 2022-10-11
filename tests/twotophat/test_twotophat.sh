#!/bin/bash

#========================================================================
# a bash script to run MICA on test data
# 
# by Yan-Rong Li, Oct 11, 2022
#========================================================================
# create folders
if [ ! -d "data/" ]; then
  mkdir data
fi
cp sim_data.txt data/

# temporarily backup param file
if [ ! -d "param/" ]; then
  mkdir param
fi
cp param_twotophat param/

# copy mica exectuable 
cp ../../mica2 .

# run mica 
mpiexec -n 6 ./mica2 param/param_twotophat

# copy plotfig.py and plot results
cp ../../plotfig.py .
python plotfig.py --param param/param_twotophat --resp_input resp_input.txt

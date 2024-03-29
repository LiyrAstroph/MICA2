************************
Resuming from Last Run
************************

``mica2`` supports resuming from last run. The adopted DNest sampling code ``CDNest`` (see :ref:`cdnest_label`)
outputs a number of sampling snapshots with a name as ``restart_dnest.txt_xxx`` for continuum reconstruction 
and ``restart_dnest1d.txt_yy_xxx``
for time delay analysis (``mica2`` first reconstructs the continuum light curve to determine the best values of parameters 
for the damped random walk model and uses these values as priors for subsequent time delay analysis.). Here, 
`xxx` mean the steps when the snapshot is saved and `yy` mean the number of Gaussians.

The default steps for continuum reconstruction are 2000, so the corresponding snapshots name is ``restart_dnest.txt_2000``. 
The default steps for time delay analysis are 2000 too, so the corresponding snapshots name is ``restart_dnest1d.txt_yy_2000``. 
If one sets different steps, the snapshot name is changed accordingly. 

Rename the above files as 

.. code:: bash

  cp restart_dnest.txt_xxx restart_dnest.txt
  cp restart_dnest1d.txt_yy_xxx restart_dnest1d.txt_yy     # namely, remove "_xxx"

Then change the step numbers in the option files (e.g., param/OPTIONS1D) and resume from last run as 

.. code:: bash 
  
  mpiexec -n np ./mica2 para/param param/OPTIONS1D -r 

where "np" is the number of cores. It must be the same as the number of cores used in last run.

In the Python verson `pymica`, a function `restart()` is provided to do resuming from last run. See 
the Python script `example.py` in the subfolder `tests/python` in the source package.
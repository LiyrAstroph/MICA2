************************
Resuming from Last Run
************************

``mica2`` supports the functionality of resuming from last run. The adopted DNest sampling code ``CDNest`` (see :ref:`cdnest_label`)
outputs a number of sampling snapshots with a name as ``restart_dnest.txt_xxx`` for continuum reconstruction and ``restart_dnest1d.txt_xxx``
for time delay analysis (``mica2`` first reconstructs the continuum light curve to determine the best values of parameters 
for the damped random walk model and uses these values as priors for subsequent time delay analysis.)

The default steps for continuum reconstruction are 1500, so the corresponding snapshots name is ``restart_dnest.txt_1500``. 
The default steps for time delay analysis are 1500 too, so the corresponding snapshots name is ``restart_dnest1d.txt_1500``. If one sets 
different steps, change the snapshot name accordingly. 

Rename the above files as 

.. code:: bash

  cp restart_dnest.txt_xxx restart_dnest.txt
  cp restart_dnest1d.txt_xxx restart_dnest1d.txt

Then changes the step numbers in the option files (e.g., param/OPTIONS1D) and resume from last run as 

.. code:: bash 
  
  mpiexec -n np ./mica2 para/param param/OPTIONS1D
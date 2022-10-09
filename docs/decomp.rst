************************
Decomposing Light Curves
************************

``mica2`` can decompose the echoing light curve into components corresponding to each Gaussian component in the transfer function, 
via supplying an option "-d" and rerun ``mica`` based on previously obtained posterior samples, 

.. code:: bash

   mpiexec -n np ./mica2 param/param        # normal mica analysis and generate posterior samples
   mpiexec -n np ./mica2 param/param -d     # decompose the light curve

Note that the second command does not redo sampling, instead, it only decompose the light curve.
The results are output to files named as "data/pline.txt_xx_comp1", where "xx" is the number of Gaussians and "comp1" means 
the first Gaussian component. 

``mica2`` provides a python script to plot the decomposed light curves. Supply the script the option file (param/param) and run 
the script as

.. code:: bash 

   python plotdecomp.py param/param


.. figure:: _static/fig_decomp.jpg
  :scale: 50 %
  :align: center
  
  (Top) The  transfer function with two Gaussians. (Center) Continuum light
  curve. (Bottom) Emission line light curve and its decomposition into two 
  components according to the two-component transfer function.
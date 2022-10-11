.. _tests_label:

*****
Tests
*****

Make sure that your MICA compiling is successful. After that, change to the subdirectory ``tests``.
There are two tests provided: ``twogauss`` and ``twotophat``. As their names indicate, they are for 
two-Gaussian and two-top-hat transer functions, respectively. Change to any one folder and execute
the bash script 

.. code:: bash 

    ./test_twogauss.sh 

or 

.. code:: bash

    ./test_twotophat.sh

The script will create folders ``data`` and ``param``, place appropriate files in the two folders, copy 
the executable file ``mica2`` and ploting script ``plotfig.py``, and then run ``mica2`` and plot the results.

Here is the output for the twogauss tests.

.. figure:: _static/fig_twogauss.jpg
  :scale: 30 %
  :align: center

  Results using two Gaussians. (Left) transfer functions; (Right) Fits to light curves.

The obtained evidences are::

    # number_of_gaussians     ln(evidence)
    1       673.680528
    2       814.866840

It is clear that the case of two Gaussians are decisively perferred over the case of one Gaussian. 
The obtained transfer function is remarkably consistent with the input.

.. note::

    MICA does not take into account **a constant normalization factor** for the likelihood function,
    so the obtained Bayesian evidence might be larger than 1. In this case, only the differences in 
    evidence make sense.
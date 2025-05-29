.. sys_err_label:

************************
Systematic Error Parameter
************************

In some cases the data errors are underestimated. MICA can add a systematic error parameter 
to augment the data errors. MICA parameterizes the systematic error as 

.. math::
  x = \ln\left(1+\frac{\sigma_{\rm sys}}{\sigma_{\rm data, mean}}\right),

where :math:`\sigma_{\rm sys}` is the systematic error and :math:`\sigma_{\rm data, mean}` is 
the mean of the data errors, defined as 

.. math:: 
  \sigma_{\rm data, mean} = \frac{1}{N}\sum_{i=1}^N \sigma_i,

where :math:`\sigma_i` is the data error of the :math:`i^{\rm th}` data point and :math:`N` is the 
number of data points.

MICA outputs the parameters list in the file **data/para_names_line.txt_xx**, where **xx**
represents the number of components used. This file looks like::

  0 sys_err_con               LOGUNI   0.000000   0.000000    1    0.000000e+00
  1 sigmad                    LOG     -4.184615  -3.473817    0  -1.797693e+308
  2 taud                      LOG      0.000000   7.359338    0  -1.797693e+308
  3 sys_err_line              LOGUNI   0.000000   2.397895    1    0.000000e+00
  4 0-th_component_amplitude  LOG    -11.512925   2.302585    0  -1.797693e+308
  5 0-th_component_center     UNI      0.000000 100.000000    0  -1.797693e+308
  6 0-th_component_sigma      LOG     -1.609438   4.422849    0  -1.797693e+308

As the name indicates, the 0th and 3rd parameters are the systematic error parameters 
for the continuum and line datasets, respectively. Then one gets the posterior samples 
of these parameters from the file **data/posterior_sample1d.txt_xx**. For example, with 
Python, one can do as follows

.. code-block:: python

    import numpy as np
    sample = np.loadtxt('data/posterior_sample1d.txt_1')
    
    err_con_mean = np.mean(con_data_err[:]) # Mean of data errors for continuum
    err_line_mean = np.mean(line_data_err[:]) # Mean of data errors for line

    sys_err_con = (np.exp(sample[:, 0])-1.0) * err_con_mean # Systematic error of continuum
    sys_err_line = (np.exp(sample[:, 3])-1.0) * err_line_mean # Systematic error of line
    
    print('Systematic error of continuum:', sys_err_con)
    print('Systematic error of line:', sys_err_line)

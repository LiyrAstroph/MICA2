**************************
Frequently Asked Questions
**************************

- **How to speed up MICA running?**

  MICA invovles matrix manipulations in optimizing the posterior probablity. The time complexity scales 
  with the number of data points as :math:`O(N^3)`, so a large number of data points will cause slow running.

  One way to reudce running time is rebinning the light curve. Provided the rebinning interval is smaller than 
  the typical time lag, this will not affect the result.

- **How to determine the best number of components?**
  
  MICA outputs Bayesian evidence (in natural log base) for each number of components into the file ``evidence.txt``. 
  One can use the Bayesian factor to select the best number of components. Simply speaking, the case with the largest
  evidence is preferred.

- **How is the systematic error parameter defined?**

  MICA parameterizes the systematic error as 

  .. math::
    x = \ln\left(1+\frac{\sigma_{\rm sys}}{\sigma_{\rm data, mean}}\right),

  where :math:`\sigma_{\rm sys}` is the systematic error and :math:`\sigma_{\rm data, mean}` is 
  the mean of the data errors, defined as 

  .. math:: 
    \sigma_{\rm data, mean} = \frac{1}{N}\sum_{i=1}^N \sigma_i,

  where :math:`\sigma_i` is the data error of the :math:`i^{\rm th}` data point and :math:`N` is the 
  number of data points. See :ref:`sys_err_label` for more details.

- **How to change the time ranges of light curve reconstruction?**
  
  The time ranges can be modified by the options `TimeRecLowExt` and `TimeRecUppExt`. For example, 
  if denote the origin time rage as [t1, t2], the options will change the time range into 
  [t1+TimeRecLowExt, t2+TimeRecUppExt].
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
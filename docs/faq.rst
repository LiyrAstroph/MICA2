**************************
Frequently Asked Questions
**************************

- **How to speed up MICA running?**

  MICA invovles matrix manipulations in optimizing the posterior probablity. The time complexity scales 
  with the number of data points as :math:`O(N^3)`, so a large number of data points will cause slow running.

  One way to reudce running time is rebinning the light curve. Provided the rebinning interval is smaller than 
  the typical time lag, this will not affect the result.
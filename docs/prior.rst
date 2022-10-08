********************
Prior Types of Lags
********************

``mica2`` provides four types of lag (namely, Gaussian center) priors by specifying the option **TypeLagPrior** in the parameter file.

* **Type 0**
  
  The lags overall have a prior range of :math:`(T_0, T_1)`, which are specified by **LagLimitLow** and **LagLimitUpp** options. 
  They obey the restrict that :math:`T_0 < \tau_0 < \tau_1 < \tau_2 <... < T_1`.

* **Type 1**

  The lags have prior ranges as:

  :math:`T_0 + 0*W < \tau_0 < T_0 + 1*W`

  :math:`T_0 + 1*W < \tau_1 < T_0 + 2*W`

  ...
  
  :math:`W = (T_1 - T_0)/K`, where :math:`K` is the number of Gaussians.

* **Type 2**
  
  The lags are fixed at specific values and there is no limit on Guassian sigma.

  :math:`\tau_0 = T_0 + 0*\Delta T`

  :math:`\tau_1 = T_0 + 1*\Delta T`

  ...

  :math:`\Delta T = (T_1 - T_0)/(K-1)`, where :math:`K` is the number of Gaussians.

* **Type 3**

  The lags are fixed at specific values and Gaussian sigma ranges at :math:`(\Delta T/2, \Delta T)`.
  
  :math:`\tau_0 = T_0 + 0*\Delta T`

  :math:`\tau_1 = T_0 + 1*\Delta T`

  ...

  :math:`\Delta T = (T_1 - T_0)/(K-1)`, where :math:`K` is the number of Gaussians.

For Types 2 and 3, it is better to use a relatively large number of Gaussians.

**Note that for top-hat transfer function, the above types are generally similar 
(amplitude, center, and sigma of Gaussians correspond to amplitude, center, and width of top-hats, 
respectively), except for type 3, in which 
the top-hat widths are fixed to be** :math:`\Delta T/2`.
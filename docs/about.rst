
*********************
What is MICA?
*********************
``mica2`` is a non-parameteric approach to analyze light curves in reverberation mapping and infer the transfer functions. 
A transfer function or delay map relates a time series  to its driving time series as

.. math::
  
  L(t) = \int \Psi(\tau) C(t-\tau) d\tau.

``mica2`` expresses the transfer function into a family of displaced Gaussians,

.. math::

  \Psi(\tau) = \sum_{k=1}^{K} \frac{f_k}{\sqrt{2\pi}\omega_k} \exp\left[-\frac{(\tau-\tau_k)^2}{2\omega_k^2}\right].

.. note::
  There is a factor :math:`1/\sqrt{2\pi}\omega_k` before the exponential 
  in the above transfer function.

.. figure:: _static/fig_sch_loc.jpg
  :scale: 50 %
  :align: center
  
  Schematic of the transfer function for a system that consists of discrete clouds.

As an alternative option, ``mica2`` also supports top-hat transfer functions as 

.. math::

  \Psi(\tau) = \sum_{k=1}^{K} \frac{f_k}{2\omega_k} H(\tau, \tau_k, \omega_k),

where :math:`H(\tau, \tau_k, \omega_k)` is the top-hat function

.. math:: 

  H(\tau, \tau_k, \omega_k) =~1~{if}~\tau_k-\omega_k \leqslant \tau \leqslant \tau_k + \omega_k

                            =~0~else~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an example for reverberation mapping analysis of the light curves from Hu et al. (2020) using **MICA**,

.. figure:: _static/fig_pg2130.jpg
  :scale: 30 %
  :align: center

  Reverberation mapping analysis of the light curve data for PG 2130+099 (Hu et al. 2020, ApJ, 890, 71).
  The left panel shows transfer functions and the right panel shows the light curves and their reconstructions.


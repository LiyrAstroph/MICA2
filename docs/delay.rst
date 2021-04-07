********************
Measuring Time Lags
********************
**MICA** outputs the posterior sample of parameters into a file with a name as 
**data/posterior_sample1d.txt_xx**, where **xx** represents the number of Gaussians used.
In this file, the order of parameters  is arranged as: 

(systematic error of continuum, sigmad, taud) * number of datasets

(systematic error of line, (Gaussian amplitude, center, sigma) * number of Gaussians * number of line datasets) * number of datasets


The Gaussian centers can be regarded as the time lags of the corresponding Gaussians. 
If only using one Gaussian, the Gaussian center can be directly used as the time lag. 
If using multiple Gaussians, one may take the Gaussian amplitude into account to 
compute a weighted time lag, e.g., 

.. math::
  
  \tau &= \int \Psi(\tau) \tau d\tau \\
       &= \sum_k \int \frac{f_k }{\sqrt{2\pi}\omega_k} \exp\left[-\frac{(\tau-\tau_k)^2}{2\omega_k^2}\right] \tau d\tau\\
       &= \sum_k f_k \tau_k
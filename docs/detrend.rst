************************
Detrending Light Curves
************************
In some cases, there may appear non-echoed long-term trends in light curves. ``mica2`` can account for those non-echoed 
trends by using a polynomial to do detrending. The order of the polynomial is specified in the parameter file via::
  
  FlagLongtermTrend         0                # Longterm trend in light curves, use a polynomial to fit 
                                             # input the order of the polynomial, e.g.,
                                             # 0, constant  (default)
                                             # 1, linear line 
                                             # 2, conic line
                                             # Use the default if you do not know this.


The default is using a constant and if there are visible non-echoed trends, usually a line polynomial is sufficient.

Below see an example for detrending light curves of 3C 273 from Zhang et al. (2019, ApJ, 876, 49).

.. figure:: _static/fig_3c273_detrend.jpg
  :scale: 50 %
  :align: center

  Time lag analysis of light curves of 3C 273 with long-term detrending (``FlagLongtermTrend=1``). Dashed 
  lines represent the long-term trends. 

*************************
Diffusive Nested Sampling
*************************

``mica2`` employ the diffusive nested sampling technique  to explore the posterior probability distribution, 
which was developed by Brewer et al. (2009). We wrote a C version of the code DNest developed by Brewer et al.
and made some modification for our purpose.

To initiate the sampling, we need to input options that control the sampling configurations. The options 
should be written in a text file： **OPTIONS1D**. This option file looks like::

  # File containing parameters for DNest
  # Put comments at the top, or at the end of the line.
  # Do not change the order of lines.
  # Lines beginning with '#' are regarded as comments.
  
  NumberParticles           2	     # Number of particles
  NewLevelIntervalFactor    2	     # new level interval
  ThreadStepsFactor         2     # threadSteps - how many steps each thread should do independently before communication
  MaxNumberSaves            10000	 # Maximum number of saves
  PTol                      0.1    # Likelihood tolerance in loge
  
  # Full options are:
  # NumberParticles          2
  # NewLevelIntervalFactor   2
  # SaveIntervalFactor       2
  # ThreadStepsFactor        2
  # MaxNumberLevels          0
  # BacktrackingLength       10.0
  # StrengthEqualPush        100.0
  # MaxNumberSaves           10000
  # PTol                     0.1


In generic, the most import options are the likelihood tolerance in loge (PTol) and the maximum number of saves 
(MaxNumberSaves). Unfornatuately, there is not yet a satisfactory rule to determine the best values for these options 
before running the code. If one finds the results not good or the effective sample two few, increase the maximum number of saves.
If this does not work, then decrease PTol option.

After MICA2 running finishes, one can run the Python script "postprocess.py" provided in the package as

.. code:: python

  python postprocess.py param

This will generate a PDF file **dnest_xx.pdf** in subdirectory **data/** that can be used to diagnose the 
whether option settings are appropriate. The file shows plots like 

.. figure:: _static/fig_dnest_post.jpg
  :scale: 30 %
  :align: center

  (Top) the likelihood (:math:`\log L`) and (bottom) posterior weigth (:math:`L\Delta X`) as 
  a function of prior mass :math:`X`. **For a good run, there should appear a peak in the 
  bottom panel** (see also the documentation of DNest3, https://github.com/eggplantbren/DNest3).

.. figure:: _static/fig_post_bad.jpg
  :scale: 60%
  :align: center 

  **An example for a bad run: there is no peak in the bottom panel.**

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
  
  2	    # Number of particles
  100	# new level interval
  100	# save interval
  20	  # threadSteps - how many steps each thread should do independently before communication
  60	  # maximum number of levels
  10	  # Backtracking scale length (lambda in the paper)
  100	  # Strength of effect to force histogram to equal push. 0-10 is best. (beta in the paper)
  1500	  # Maximum number of saves


In generic, the most import options are the maximum number of levels and the maximum number of saves. Unfornatuately, 
there is not yet a satisfactory rule to determine the best values for these options before running the code. If one 
finds the results not good or the effective sample two few, increase the maximum number of saves. If this does not 
work, then increase the maximum number of levels.
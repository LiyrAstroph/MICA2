****************
Seasonal Gap
****************

For multi-season campaigns, the unavoidable seasonal gaps might bias the reverberation mapping analysis.
``mica2`` can cope with seasonal gaps by excluding the regions in parameter space that might be connected 
with gaps. Those regions are defined as 

.. math:: 
  
  |\tau_k - \tau_{\rm gap} - n\times 365.25| < \omega_{\rm gap},~~~~ \omega_k < \omega_{\rm gap},

where :math:`n` is any integer,  :math:`\tau_k` and :math:`\omega_k` are the time lag and width of :math:`k`-th component, :math:`\tau_{\rm gap}`
and :math:`\omega_{\rm gap}` are the central time lag resulting from gaps and gap width, respectively. By default, 
``mica2`` sets 

.. math:: 

  \tau_{\rm gap} = 365.25/2,

and :math:`\omega_{\rm gap}` is the assigned as the mean gap width of the data.

To turn on this functionality, edit the option in the parameter file::

  FlagGap                 0                  # whether include seasonal gap
                                             # 0: no; 1: yes.
                                             # default: 0

If the default values are not satifactory, edit the option::

  #StrGapPrior            [182.6:140.0]      # gap priors if the default priors are not good enough
                                             # valid when FlagGap == 1
                                             # format: [gap_center_set1:gap_width_set1:gap_center_set2:gap_width_set2...]
                                             # gap_center_set1: gap center for 1st dataset (+n*year will also be included)
                                             # gap_width_set1:  gap width for 1st dataset
                                             # default: None

In the Python version, use the arguments as 

.. code-block:: python
  
  model = pymica.gmodel()
  model.setup(data=data_input, ..., flag_gap=True)

or input the desired gap information as

.. code-block:: python
  
  model = pymica.gmodel()
  model.setup(data=data_input, ..., flag_gap=True, gap_prior=[[182.625, 100],])

where ``gap_prior`` is a list and specifies the central time lag (182.625 day) and width (100 day) of gaps for all datasets. 
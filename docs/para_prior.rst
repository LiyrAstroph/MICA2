.. _para_prior_label:

*************************************
Adjusting Prior Ranges of Parameters
*************************************

MICA sets default prior ranges for parameters in MCMC sampling. In some cases, these default prior 
ranges might be inapproprite.  MICA allows to adjust prior ranges as follows.

First, generate the default prior ranges adopted in MICA with the command (alternatively, one can also use 
the generated prior file from the last run)

.. code-block:: bash

    ./mica2 param -n 


This will output the prior ranges to a file named `data/para_names_xx.txt`. The content looks like::

    0 sys_err_con               LOGUNI   0.000000   0.000000    1    0.000000e+00
    1 sigmad                    LOG     -4.160532  -3.498325    0  -1.797693e+308
    2 taud                      LOG      1.000000   7.359338    0  -1.797693e+308
    3 sys_err_line              LOGUNI   0.000000   2.397895    1    0.000000e+00
    4 0-th_component_amplitude  LOG    -11.512925   2.302585    0  -1.797693e+308
    5 0-th_component_center     UNI      0.000000 100.000000    0  -1.797693e+308
    6 0-th_component_sigma      LOG     -1.609438   4.422849    0  -1.797693e+308


Here, the columns are respectively the ID, Name, Prior type, Min, Max, Fix, and Val.
`Min` and `Max` columns represent the lower and upper limits. For `Fix` column, `0` means not fixed and `1` means fixed. 
When the parameter is fixed, the `Val` columns represent the fixed value. The lines starting with "#" will be neglected.

**Note that do not change the format the file, otherwise, there will be an error when reading in it.**

If the prior range of some parameters needs to change, just adjust the corresponding `Min` and `Max` columns. 
Then save the edited file to a new file, e.g., say, `data/new_prior.txt`. 
Note that for the sake of brevity, one can only keep those lines for the parameters to be adjusted. The rest lines can be removed. 
However, it is still fine to keep all lines. For example, one can edit the 1st and 3rd parameters (counting from 0) as::

    1 sigmad                    LOG     -4.160532  -3.498325    0  -1.797693e+308
    2 taud                      LOG      1.000000   7.359338    0  -1.797693e+308

Afterwards, pass this new prior file to MICA as 

.. code-block:: bash

    mpiexec -n 6 ./mica2 param -l data/new_prior.txt

MICA will read in the prior ranges and use them for MCMC sampling.


In the Python version, one can generate the default prior ranges as 

.. code-block:: python

    model = pymica.gmodel()
    ...
    model.print_para_names()

To load priors from a file, say, `data/new_prior.txt`, call the function 

.. code-block:: python 

    model.set_priors("data/new_prior.txt")
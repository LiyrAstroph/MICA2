.. _format:

***************************
MICA-Format Input Data
***************************

MICA accepts input data in a specific format as described in :ref:`getting_started`. 
The Python version of MICA ``pymica`` provide a function ``format_data()`` to convert 
the input data into the required format.

If one has a driving light curve and a list of responding light curves, call the function as follows.

.. code-block:: python

    from pymica.utility import format_mica

    # generate an output file "filename_output.txt".
    # each light curve "lc" is a numpy array, with a size of n*3, namely, containing 3 columns: 
    # time, flux, and flux error.

    data1 = lc_driving
    data2 = [lc1, lc2, lc3]
    format_mica("filename_output.txt", data1, data2)   

If one has several datasets with each containing a driving light curve and a list of responding 
light curves, call the function as follows.

.. code-block:: python

    from pymica.utility import format_mica

    # generate an output file "filename_output.txt".
    # each light curve "lc" is a numpy array, with a size of n*3, namely, containing 3 columns: 
    # time, flux, and flux error.

    data1 = [lc_driving1, lc_driving2, lc_driving3]
    data2 = [[lc11, lc12, lc13], [lc21, lc22], [lc31, lc32, lc33]]
    format_mica("filename_output.txt", data1, data2)   

In the ``vmap`` mode, there is no driving light curve. One can call the function as follows.

.. code-block:: python

    from pymica.utility import format_mica

    # generate an output file "filename_output.txt".
    # each light curve "lc" is a numpy array, with a size of n*3, namely, containing 3 columns: 
    # time, flux, and flux error.

    data1 = None
    data2 = [lc1, lc2, lc3]
    format_mica("filename_output.txt", data1, data2)   


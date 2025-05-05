***************************************************
Transfer Function and Its Uncertainties 
***************************************************

``mica2`` calculates the transfer function and its uncertainties from 
the posterior samples of parameters as follows. Given one set of parameters, 
one can calculate the corresponding transfer function. Looping over 
the whole posterior sample with :math:`N` sets of parameters, one get :math:`N` 
transfer functions. Then one can assign the median as the best transfer function,
and use the 68.3% quantiles to determine the uncertainties.

Here is the Python procedure used in ``mica2``.

.. code-block:: python
    
    sample = np.atleast_2d(np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau))

    if typetf == 0: # gaussian
        for i in range(sample.shape[0]):
          # loop over gaussians
          if typemodel == 0:  # general model
            for k in range(ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)

          elif typemodel == 1: # pmap model 
            k = 0
            amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
            tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)
            for k in range(1, ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)

    else:  # tophats
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0: # general model
            for k in range(ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))

          elif typemodel == 1: # pmap model 
            k = 0
            amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
            tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))
            for k in range(1, ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))
      
    tran_best = np.percentile(tran, 50.0, axis=0)
    tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
    tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)
    
    ax.plot(tau, tran_best, color='k')
    ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')
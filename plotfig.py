# MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
#
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os
import configparser as cp 
import argparse

__all__ = ["plot_results"]

def plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend, flagnegresp, typetf, typemodel, resp_input, 
                 doshow=True, tf_lag_range=None, hist_lag_range=None):
  """
  reconstruct line lcs according to the time sapns of the continuum.
  """
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif', size=15)

  sample = np.atleast_2d(np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau))
  data = np.loadtxt(fdir+fname)
  sall = np.loadtxt(fdir+"/data/pall.txt_%d"%ngau)

  if flagtrend > 0:
    trend = np.loadtxt(fdir+"/data/trend.txt_%d"%ngau)

  fp = open(fdir+fname)
  # read numbe of datasets
  line = fp.readline()
  nd = int(line[1:])
  if flagvar == 1:
    num_params_var = 3
  else:
    num_params_var = 3*nd
  
  # number of parameters for long-term trend
  nq = flagtrend + 1
  
  # read number of data points in each dataset
  nl = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = np.array([int(i) for i in ls])
    nl.append(ns)
  fp.close()

  # read number of points of reconstructions
  fp = open(fdir+"/data/pall.txt_%d"%ngau, "r")
  line = fp.readline()
  nl_rec = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = np.array([int(i) for i in ls])
    nl_rec.append(ns)
  fp.close()

  # assign index of cont data
  indx_con_data = []
  indx_con_rec = []
  indx_con_data.append(0)
  indx_con_rec.append(0)
  for i in range(1, nd):
    ns = nl[i-1]
    ns_rec = nl_rec[i-1]
    indx_con_data.append(np.sum(ns) + indx_con_data[i-1])
    indx_con_rec.append(np.sum(ns_rec) + indx_con_rec[i-1])

  # assign index of the parmaeter for the first line of each dataset 
  indx_line = []
  indx_line.append(num_params_var)
  for i in range(1, nd):
    if flagtran == 1:
      indx_line.append(num_params_var)
    else:
      indx_line.append(indx_line[i-1] + (len(nl[i-1])-1)*(1+ngau*3))

  # print time lags, median, and 68.3% confidence limits
  if flagnegresp == False:
    print("========No. Gaussian/Tophat: %d========="%ngau)
    print("ID:      lag     -elo      +eup")
    sample_lag = np.zeros(sample.shape[0])
    weight_lag = np.zeros(sample.shape[0])
    #loop over datasets
    for m in range(nd):
      print("Dataset %d"%m)
      #loop over lines
      ns = nl[m]
      for j in range(1, len(ns)):
        sample_lag[:] = 0.0
        weight_lag[:] = 0.0
        if typemodel == 0 or typemodel == 2:  # general, vmap model
          for k in range(ngau):
            if flagnegresp == 0:  # no negative response
              sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] * np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              weight_lag[:] +=  np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            else:
              sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              weight_lag[:] +=  1.0

        elif typemodel == 1: # pmap model
          k = 0
          sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] * np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          weight_lag[:] +=  np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          for k in range(1, ngau):
            sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] * np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            weight_lag[:] +=  np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                                    sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
    
        lag, err1, err2 = np.quantile(sample_lag/weight_lag, q=(0.5, (1.0-0.683)/2.0, 1.0-(1.0-0.683)/2.0))
        err1 = lag-err1
        err2 = err2 - lag
        print("Line %d: %.3f -%.3f +%.3f"%(j, lag, err1, err2))
  
  # load input respon function
  if resp_input != None:
    tran_input = np.loadtxt(resp_input)

  dtau = tau_upp - tau_low 
  ntau = 500
  tran = np.zeros((sample.shape[0], ntau))
  
  shift = 0.0
  
  # open pdf file
  pdf = PdfPages(fdir+"/data/fig_%d.pdf"%ngau)
  print("Plotting to %s."%(fdir+"/data/fig_%d.pdf"%ngau))
  
  # set x-axis coordinates of sub figures
  if flagnegresp == False:
    figlc_center = 0.05
    figlc_centroid = 0.22
  else:
    figlc_center = 0.17
    # figlc_centroid = 0.22

  idx_q = 0 # index for long-term trend parameters
  for m in range(nd):
    ns = nl[m]
    ns_rec = nl_rec[m]
    fig = plt.figure(figsize=(12, 4+2*(len(ns)-2)))
    
    #===================================
    # plot continuum
    con0 = data[indx_con_data[m]:indx_con_data[m]+ns[0], :]
    sall_con0 = sall[indx_con_rec[m]:(indx_con_rec[m]+ns_rec[0]), :] 
    
    axheight = 0.8/(len(ns))
    ax = fig.add_axes((0.56, 0.95-axheight, 0.35, axheight))
    
    ax.errorbar(con0[:, 0]-shift, con0[:, 1], yerr=con0[:, 2], ls='none', color='b', zorder=10, marker='o', markersize=1.5, elinewidth=0.5)
    
    ax.plot(sall_con0[:, 0]-shift, sall_con0[:, 1], color='k')
    ax.fill_between(sall_con0[:, 0]-shift, y1=sall_con0[:, 1]-sall_con0[:, 2], y2=sall_con0[:, 1]+sall_con0[:, 2], color='darkgrey')
    
    # plot long-term trend
    if flagtrend > 0:
      xlim = ax.get_xlim()
      x = np.linspace(sall_con0[0, 0], sall_con0[-1, 0], 100)
      y = np.zeros(100)
      for j in range(nq):
        y += trend[idx_q + j, 0] * x**(j)
  
      ax.plot(x, y, ls='--', color='grey')

    ax.xaxis.set_tick_params(labeltop=False)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    ax.yaxis.set_tick_params(labelright=True)
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('Flux')
    ax.minorticks_on()
    #ax.xaxis.set_label_position("top")
    #ax.set_xlabel('HJD - 2450000')

    # set ylim
    # note in vmap case, no continuum data
    if con0.shape[0] > 0:
      ymin = np.min(con0[:, 1])
      ymax = np.max(con0[:, 1])
      dy = ymax - ymin
      ax.set_ylim(ymin-0.1*dy, ymax+0.1*dy)
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xlim0 = xlim # restore the continuum time range
  
    # plot line
    # set time lag range for Gaussian centers and centriods
    if hist_lag_range is None:
      tau1 = 1.0e10
      tau2 = -1.0e10
      tau1_cent = 1.0e10
      tau2_cent =-1.0e10
      if typetf in [0, 1]: 
        for j in range(1, len(ns)):    
          for k in range(ngau):
            tau1 = np.min((tau1, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
            tau2 = np.max((tau2, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
        
        tau1_cent = tau1
        tau2_cent = tau2 
      elif typetf == 2:
        for j in range(1, len(ns)):  
          for k in range(ngau): # gamma use peak
            tau1 = np.min((tau1, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                                      +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
            tau2 = np.max((tau2, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                      +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))

            tau1_cent = np.min((tau1_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
              + 2.0*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
            tau2_cent = np.max((tau2_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
              + 2.0*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))
      elif typetf == 3: 
        for j in range(1, len(ns)):  
          for k in range(ngau): # exp use peak
            tau1 = np.min((tau1, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1], q=0.0)))
            tau2 = np.max((tau2, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1], q=1.0)))

            tau1_cent = np.min((tau1_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
              + np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
            tau2_cent = np.max((tau2_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
              + np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))
    else:
      tau1 = hist_lag_range[0]
      tau2 = hist_lag_range[1]

      tau1_cent = hist_lag_range[0]
      tau2_cent = hist_lag_range[1]


    # set time lag range for transfer function
    if tf_lag_range is None: 
      tau1_tf = 1.0e10
      tau2_tf = -1.0e10
      for j in range(1, len(ns)): 
        if typetf == 0:  # gaussian  
          for k in range(ngau):
            tau1_tf = np.min((tau1_tf, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              -3*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
            tau2_tf = np.max((tau2_tf, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              +3*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
        elif typetf == 1: # tophats
          for k in range(ngau):
            tau1_tf = np.min((tau1_tf, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              -1.5*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
            tau2_tf = np.max((tau2_tf, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              +1.5*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))))
        elif typetf == 2:  # gamma
          for k in range(ngau):
            tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              -0.2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
            tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              +6*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))
        else:  # exp
          for k in range(ngau):
            tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              -0.2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
            tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                              +6*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))

      tau1_tf = np.min((tau_low, tau1_tf))
      tau2_tf = np.max((tau_upp, tau2_tf))
    else:
      tau1_tf = float(tf_lag_range[0])
      tau2_tf = float(tf_lag_range[1])
        
    # now do plotting
    for j in range(1, len(ns)):
      hb = data[indx_con_data[m] + np.sum(ns[:j]):indx_con_data[m] + np.sum(ns[:j+1]), :] 
      sall_hb = sall[(indx_con_rec[m] + np.sum(ns_rec[:j])):(indx_con_rec[m] + np.sum(ns_rec[:j+1])), :]
      
      # histogram of time lags 
      ax = fig.add_axes((figlc_center, 0.95-(j+1)*axheight, 0.16, axheight))

      for k in range(ngau):
        if typetf in [0, 1]: # gaussian or tophat
          cen = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        elif typetf == 2:  # gamma, use peaks
          cen =  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        elif typetf == 3:  # exp, use peaks
          cen =  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]

        cen_min = np.min(cen)
        cen_max = np.max(cen)
        bins = np.max((5, int((cen_max-cen_min)/(tau2-tau1 + 1.0e-100) * 20)))
        if k == 0:
          ax.hist(cen, density=True, bins=bins, alpha=1)
        else:
          ax.hist(cen, density=True, bins=bins, alpha=0.6)

      ax.set_xlim((tau1-(tau2-tau1)*0.1, tau2+(tau2-tau1)*0.1))
      
      # vmap model, no need to plot for the first lc, which has a zero lag wrt the driving lc. 
      if typemodel == 2 and j == 1:
        ax.set_visible(False)
      
      if flagnegresp == False:
        ax.yaxis.set_tick_params(labelleft=False)

      if (typemodel != 2 and j == 1) or (typemodel == 2 and j == 2):
        if typetf in [0, 1]:
          ax.set_title("Centers")
        else:
          ax.set_title("Peaks")

      ax.minorticks_on()
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")
      
      # only plot centroid lag when no negative response
      if flagnegresp == False:
        ax = fig.add_axes((figlc_centroid, 0.95-(j+1)*axheight, 0.16, axheight))
        
        if typemodel == 0 or typemodel == 2: # centroid time lag
          cent = np.zeros(sample.shape[0])
          norm = np.zeros(sample.shape[0])
          for k in range(ngau):

            if flagnegresp == 0: # no negative responses
              if typetf in [0, 1]:
                norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                        * sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              elif typetf == 2:
                norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                      * (sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                        +2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
              elif typetf == 3:
                norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
                cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                      * (sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                        +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
            else:
              norm += 1.0
              cent += sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          
          ax.hist(cent/norm, density=True, range=(tau1_cent, tau2_cent), bins=20)
          ax.set_xlim((tau1_cent-0.1*(tau2_cent-tau1_cent), tau2_cent+0.1*(tau2_cent-tau1_cent)))
          ax.minorticks_on()
          ax.yaxis.set_tick_params(labelleft=False)

          if (typemodel != 2 and j == 1) or (typemodel == 2 and j == 2):
            ax.set_title("Centroid")
            
          if j != len(ns)-1:
            ax.xaxis.set_tick_params(labelbottom=False)
          else:
            ax.set_xlabel("Time Lag (day)")

        elif typemodel == 1: # ratio hist
          for k in range(1, ngau):
            ratio = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3]/np.log(10.0)
            if k == 1:
              ax.hist(ratio, density=True, bins=20, alpha=1)
            else:
              ax.hist(ratio, density=True, bins=20, alpha=0.6)

            ax.yaxis.set_tick_params(labelleft=False)
            ax.minorticks_on()

            if j != len(ns)-1:
              ax.xaxis.set_tick_params(labelbottom=False)
            else:
              ax.set_xlabel("$\log R$")

            if j == 1:
              ax.set_title("Response Ratio")
        
        # vmap model, no need to plot for the first lc, which has a zero lag wrt the driving lc. 
        if typemodel == 2 and j == 1:
          ax.set_visible(False)

      # transfer function
      ax = fig.add_axes((0.39, 0.95-(j+1)*axheight, 0.16, axheight))
      
      tau = np.linspace(tau1_tf, tau2_tf, ntau)
      tran[:, :] = 0.0
      if typetf == 0: # gaussian
        for i in range(sample.shape[0]):
          # loop over gaussians
          if typemodel == 0 or typemodel == 2:  # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)

          elif typemodel == 1: #pmap model 
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

      elif typetf == 1:  # tophats
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0 or typemodel == 2: # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))

          elif typemodel == 1: #pmap model 
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
      elif typetf == 2:  # gamma
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0 or typemodel == 2: # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              idx_tau = np.where(tau >= cen)[0]
              tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)

          elif typemodel == 1: #pmap model 
            k = 0
            amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
            idx_tau = np.where(tau >= cen)[0]
            tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
            for k in range(1, ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
      else:  # gamma
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0 or typemodel == 2: # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              idx_tau = np.where(tau >= cen)[0]
              tran[i, idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)

          elif typemodel == 1: #pmap model 
            k = 0
            amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
            idx_tau = np.where(tau >= cen)[0]
            tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
            for k in range(1, ngau):
              amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)
      
      tran_best = np.percentile(tran, 50.0, axis=0)
      tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
      tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)

      out = np.column_stack((tau, tran_best, tran_best-tran1, tran2-tran_best))
      np.savetxt(fdir+f"/data/tranfunc_{m}.txt_{ngau}", out, fmt="%f")
      
      ax.plot(tau, tran_best, color='k')
      ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')

      #plot input response function
      if resp_input != None:
        tau_min = np.max((tau[0], tran_input[0, 0]))
        tau_max = np.min((tau[-1], tran_input[-1, 0]))
        
        # normalize with the same tau range
        idx_tran = np.where((tau>=tau_min)&(tau<=tau_max))[0]
        idx_tran_input = np.where((tran_input[:, 0]>=tau_min)&(tran_input[:, 0]<=tau_max))[0]

        if flagnegresp == False:
          tran_scale = np.sum(tran_best[idx_tran])*(tau[1]-tau[0])  \
            /(np.sum(tran_input[idx_tran_input, 1])*(tran_input[1, 0]-tran_input[0, 0]))
        else:
          tran_scale = (np.max(tran_best[idx_tran])-np.min(tran_best[idx_tran])) \
            /(np.max(tran_input[idx_tran_input, 1])-np.min(tran_input[idx_tran_input, 1]))
        
        tran_input[:, 1] *= tran_scale
        ax.plot(tran_input[:, 0], tran_input[:, 1], label='input', lw=1)
        ax.legend()
      
      # determine the best range of time lag for gamma tf
      if tf_lag_range is None:
        if typetf == 2:
          idx_best_max = np.argmax(tran_best)
          idx_best_upp = np.where(tran_best[idx_best_max:]<tran_best[idx_best_max]*0.01)[0]
          if len(idx_best_upp) == 0:
            tau_best_upp = tau[-1]
          else:
            tau_best_upp = tau[idx_best_max + idx_best_upp[-1]]
          ax.set_xlim(tau[0], tau_best_upp)
        else:
          ax.set_xlim((tau[0], tau[-1]))
      else:
        ax.set_xlim(tau1_tf, tau2_tf)

      ylim = ax.get_ylim()
      if resp_input == None:
        ax.set_ylim(ylim[0], np.min((ylim[1], np.max(np.max(tran_best)*1.5))))
      else:
        ax.set_ylim(ylim[0], np.min((ylim[1], np.max((np.max(tran_best)*1.5, np.max(tran_input[:, 1])*1.5)))))
      
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")
      
      if j == 1:
        ax.set_title("Transfer Function")
      
      if flagnegresp == False:
        ax.yaxis.set_tick_params(labelleft=False)
      ax.minorticks_on()
      #if(tau[0]<0.0):
      #  ax.axvline(x=0.0, ls='--', color='red')
      
      # then line light curve
      ax = fig.add_axes((0.56, 0.95-(j+1)*axheight, 0.35, axheight))
      
      ax.errorbar(hb[:, 0]-shift, hb[:, 1], yerr=hb[:, 2], ls='none', zorder=10, color='b', marker='o', markersize=1.5, elinewidth=0.5)
      
      ax.plot(sall_hb[:, 0]-shift, sall_hb[:, 1], color='k')
      ax.fill_between(sall_hb[:, 0]-shift, y1=sall_hb[:, 1]-sall_hb[:, 2], y2=sall_hb[:, 1]+sall_hb[:, 2], color='darkgrey')
      
      # plot long-term trend 
      if flagtrend > 0:
        xlim = ax.get_xlim()
        x = np.linspace(sall_hb[0, 0], sall_hb[-1, 0], 100)
        y = np.zeros(100)
        for k in range(nq):
          y+= trend[idx_q + j*nq + k, 0]* x**(k)
        
        ax.plot(x, y, ls='--', color='grey')

      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time")
      ax.yaxis.set_tick_params(labelleft=False)
      ax.yaxis.set_tick_params(labelright=True)
      ax.yaxis.set_label_position("right")
      ax.set_ylabel("Flux")
      ax.set_xlim(xlim0[0], xlim0[1])
      
      # set ylim
      ymin = np.min(hb[:, 1])
      ymax = np.max(hb[:, 1])
      dy = ymax - ymin
      ax.set_ylim(ymin-0.1*dy, ymax+0.1*dy)

      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      ax.minorticks_on()
    
    pdf.savefig(fig)
    
    if doshow:
      plt.show()
    else:
      plt.close()

    idx_q += len(ns) * nq
  
  pdf.close()
  return

def plot_results_all(args, param, doshow=True, tf_lag_range=None, hist_lag_range=None):
  try:
    fdir = param["FileDir"]+"/"
  except:
    raise IOError("FileDir is not set!")

  try:
    flagvar = int(param["FlagUniformVarParams"])
  except:
    flagvar = 0
  
  try:
    flagtran = int(param["FlagUniformTranFuns"])
  except:
    flagtran = 0

  try:
    flagtrend = int(param["FlagLongtermTrend"])
  except:
    flagtrend = 0

  try:  
    ngau_low = int(param["NumCompLow"])
  except:
    raise IOError("NumCompLow is not set!")

  try: 
    ngau_upp = int(param["NumCompUpp"])
  except:
    raise IOError("NumCompUpp is not set!")
  
  if param["TypeModel"] == 0:
    try:
      tau_low = float(param["LagLimitLow"])
    except:
      raise IOError("LagLimitLow is not set!")
    
    try:
      tau_upp = float(param["LagLimitUpp"])
    except:
      raise IOError("LagLimitUpp is not set!")
  else:
    tau_low = 0.0
    tau_upp = 0.0
  
  try:
    fname = param["DataFile"]
  except:
    raise IOError("DataFile is not set!")
  
  try:
    typetf = int(param["TypeTF"])
  except:
    typetf = 0
  
  try:
    typemodel = int(param["TypeModel"])
  except:
    typemodel = 0
  
  try:
    flagnegresp = int(param["FlagNegativeResp"])
  except:
    flagnegresp = 0

  for ngau in range(ngau_low, ngau_upp+1):
    plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend, flagnegresp, typetf, typemodel, args.resp_input, 
                 doshow=doshow, tf_lag_range=args.tf_lag_range, hist_lag_range=args.hist_lag_range)

def _param_parser(fname):
  """
  parse parameter file
  """
  config = cp.RawConfigParser(delimiters=' ', comment_prefixes='#', inline_comment_prefixes='#', 
  default_section=cp.DEFAULTSECT, empty_lines_in_values=False)
  
  with open(fname) as f:
    file_content = '[dump]\n' + f.read()

  config.read_string(file_content)
    
  return config['dump']

if __name__ == "__main__":
  #
  parser = argparse.ArgumentParser(usage="python plotfig.py [options]")
  parser.add_argument('--param', type=str, help="parameter file")
  parser.add_argument('--resp_input', type=str, help="str, a file storing input response function")
  parser.add_argument('--tf_lag_range', type=float, nargs='+', help="time lag range for the transfer function, e.g., --tf_lag_range 0 100")
  parser.add_argument('--hist_lag_range', type=float, nargs='+', help="time lag range for the histograms, e.g., --hist_lag_range 0 100")
  args = parser.parse_args()

  if args.param == None:
    print("Please specify paramter file!")
    print("e.g., python plotfig.py --param src/param")
    sys.exit()

  fparam = args.param
  param = _param_parser(fparam)
  
  if "TypeModel" not in param:
    param["TypeModel"] = "0"

  plot_results_all(args, param)

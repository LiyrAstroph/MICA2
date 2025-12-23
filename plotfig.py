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
from enum import Enum 

__all__ = ["plot_results"]

# set the default parameters for matplotlib
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=15)
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

class Model(Enum):
  gmodel=0
  pmap = 1
  vmap = 2
  mmap = 3
  nmap = 4

def _get_hist_timelag_range(sample, indx_line, typetf, ngau, ns, m):
  """
  get the range of time lags for hist plots
  """
  tau1 = 1.0e10
  tau2 = -1.0e10
  tau1_cent = 1.0e10
  tau2_cent =-1.0e10

  for j in range(1, len(ns)):    
    for k in range(ngau):
      tt = int(typetf[k])
      if tt in [0, 1]: 
        tau1 = np.min((tau1, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
        tau2 = np.max((tau2, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
    
        tau1_cent = tau1
        tau2_cent = tau2 
      elif tt == 2: # gamma use peak
        tau1 = np.min((tau1, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                                  +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
        tau2 = np.max((tau2, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                  +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))

        tau1_cent = np.min((tau1_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
          + 2.0*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
        tau2_cent = np.max((tau2_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
          + 2.0*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))
      elif tt == 3:# exp use peak
        tau1 = np.min((tau1, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1], q=0.0)))
        tau2 = np.max((tau2, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1], q=1.0)))

        tau1_cent = np.min((tau1_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
          + np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.0)))
        tau2_cent = np.max((tau2_cent, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
          + np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=1.0)))
  
  dtau = tau2 - tau1
  tau1 -= 0.1*dtau
  tau2 += 0.1*dtau

  dtau = tau2_cent - tau1_cent
  tau1_cent -= 0.1*dtau
  tau2_cent += 0.1*dtau

  return tau1, tau2, tau1_cent, tau2_cent

def _get_ratio_range(sample, indx_line, ngau, ns, m):
  # determine ratio range 
  ratio1 = 1.0e10
  ratio2 = 1.0e-10
  for j in range(1, len(ns)): 
    for k in range(1, ngau):
      ratio1 = np.min((ratio1, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3], q=0.05)))
      ratio2 = np.max((ratio2, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3], q=0.05)))

  ratio1 /= np.log(10.0)
  ratio2 /= np.log(10.0)

  return ratio1, ratio2

def _get_tf_timelag_range(sample, indx_line, typetf, ngau, ns, m):
  # get time lag range for transfer function
  tau1_tf = 1.0e10
  tau2_tf = -1.0e10
  for j in range(1, len(ns)): 
    for k in range(ngau):
      tt = int(typetf[k])
      if tt == 0:  # gaussian  
        tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          -3*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
        tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          +3*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))
      elif tt == 1: # tophats
        tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          -1.5*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
        tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          +1.5*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))
      elif tt == 2:  # gamma
        tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          -0.2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
        tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          +6*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))
      else:  # exp
        tau1_tf = np.min((tau1_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          -0.2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.05)))
        tau2_tf = np.max((tau2_tf, np.quantile(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
                                          +6*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]), q=0.95)))

  return tau1_tf, tau2_tf

def _calculate_tran(tau, pmodel, typemodel, typetf, ngau, flagnegresp, indx_line, id, il):
  """
  calculate transfer function given a set of parameters
  """

  tran = np.zeros(len(tau))
  m = id  # dataset index
  j = il   # line index

  for k in range(len(typetf)):
    tt = int(typetf[k])
    if tt == 0: # gaussian
      # loop over gaussians
      if typemodel in [0, 2, 3]:  # general, vmap model
        if flagnegresp == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        else:
          amp =      pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

        cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        tran[:] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)/np.sqrt(2*np.pi)

      elif typemodel == 1: #pmap model 
        if k == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[:] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)/np.sqrt(2*np.pi)
        else:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
          
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[:] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)/np.sqrt(2*np.pi)

    elif tt == 1:  # tophats
      # loop over tophats
      if typemodel in [0, 2, 3]: # general, vmap model
        if flagnegresp == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        else:
          amp =      pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

        cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        
        tran[:] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))

      elif typemodel == 1: #pmap model 
        if k == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[:] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))
        else:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
          
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[:] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))

    elif tt == 2:  # gamma
      # loop over tophats
      if typemodel in [0, 2, 3]: # general, vmap model
        
        if flagnegresp == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        else:
          amp =      pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

        cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        
        idx_tau = np.where(tau >= cen)[0]
        tran[idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)

      elif typemodel == 1: #pmap model 
        if k == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          idx_tau = np.where(tau >= cen)[0]
          tran[idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
        else:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
          
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)

    else:  # exp
      # loop over tophats
      if typemodel in [0, 2, 3]: # general, vmap model
        if flagnegresp == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        else:
          amp =      pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]

        cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        
        idx_tau = np.where(tau >= cen)[0]
        tran[idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)

      elif typemodel == 1: #pmap model 
        if k == 0:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          idx_tau = np.where(tau >= cen)[0]
          tran[idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)
        else:
          amp = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                       pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+0*3+0])
          
          cen =        pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(pmodel[indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)

  return tran 

def _get_time_lag(sample, idx, tt):
  """
  get time lag from sample for the transfer function type tt
  """
  if tt in [0, 1]:
    return sample[:, idx+1] 
  elif tt == 2:
    return sample[:, idx+1] + 2*np.exp(sample[:, idx+2])
  elif tt == 3:
    return sample[:, idx+1] + 1*np.exp(sample[:, idx+2])

def plot_results_con(fdir, fname, flagtrend, doshow=True):
  """
  plot results for continuum modeling for nmap 
  """
  fp = open(fdir+fname)
  # read numbe of datasets
  line = fp.readline()
  nd = int(line[1:])

  # read number of data points in each dataset
  ncon = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = [int(i) for i in ls]
    ncon += ns
  fp.close()
  
  # reconstruction points 
  fp = open(fdir+"/data/pcon.txt")
  # read numbe of datasets
  line = fp.readline()
  nd = int(line[1:])

  # read number of data points in each dataset
  ncon_rec = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = [int(i) for i in ls]
    ncon_rec += ns
  fp.close()

  # assign index of cont data
  indx_con_data = []
  indx_con_rec = []
  indx_con_data.append(0)
  indx_con_rec.append(0)
  for i in range(1, nd):
    ns = ncon[i-1]
    ns_rec = ncon_rec[i-1]
    indx_con_data.append(ns + indx_con_data[i-1])
    indx_con_rec.append(ns_rec + indx_con_rec[i-1])

  data = np.loadtxt(fdir+fname)
  sall = np.loadtxt(fdir+"/data/pcon.txt")

  if flagtrend > 0:
    trend = np.loadtxt(fdir+"/data/trend_con.txt")
  
  # number of parameters for long-term trend
  nq = flagtrend + 1

  # open pdf file
  pdf = PdfPages(fdir+"/data/fig_con.pdf")
  print("Plotting to %s."%(fdir+"/data/fig_con.pdf"))

  shift = 0.0
  idx_q = 0
  for m in range(nd):
    ns = ncon[m]
    ns_rec = ncon_rec[m]

    fig = plt.figure(figsize=(12, 4))
    
    con0 = data[indx_con_data[m]:indx_con_data[m]+ns, :]
    sall_con0 = sall[indx_con_rec[m]:(indx_con_rec[m]+ns_rec), :] 
    
    axheight = 0.8
    ax = fig.add_axes((0.05, 0.95-axheight, 0.75, axheight))
    
    ax.errorbar(con0[:, 0]-shift, con0[:, 1], yerr=con0[:, 2], ls='none', color='b', zorder=10, marker='o', 
                markersize=3.5, elinewidth=0.5, capsize=1.0, capthick=0.5)
    
    ax.plot(sall_con0[:, 0]-shift, sall_con0[:, 1], color='k')
    ax.fill_between(sall_con0[:, 0]-shift, y1=sall_con0[:, 1]-sall_con0[:, 2], y2=sall_con0[:, 1]+sall_con0[:, 2], color='darkgrey')
    
    ax.set_ylabel('Flux')
    ax.minorticks_on()
    ax.set_xlabel('Time')

    # plot long-term trend
    if flagtrend > 0:
      xlim = ax.get_xlim()
      x = np.linspace(sall_con0[0, 0], sall_con0[-1, 0], 100)
      y = np.zeros(100)
      for j in range(nq):
        y += trend[idx_q + j] * x**(j)
  
      ax.plot(x, y, ls='--', color='grey')
    
    idx_q += nq

    #===============================================
    # residuals
    ax_res = fig.add_axes((0.82, 0.95-axheight, 0.10, axheight))
    con0_intp = np.interp(con0[:, 0], sall_con0[:, 0], sall_con0[:, 1])
    ax_res.hist((con0[:, 1]-con0_intp)/con0[:, 2], bins=20, orientation='horizontal', density=True, color='C0')
    ax_res.yaxis.set_tick_params(labelleft=False, labelright=True)
    ax_res.yaxis.set_label_position("right")
    ax_res.set_ylabel("Res./Err.")
    ax_res.set_xlabel("Hist.")
    ax_res.set_ylim(-4, 4)
    ylim_res = ax_res.get_ylim()
    ylim_res_max = max(abs(ylim_res[0]), abs(ylim_res[1]))
    ax_res.set_ylim(-ylim_res_max, ylim_res_max)
    y = np.linspace(-ylim_res_max, ylim_res_max, 100)
    x = 1.0/np.sqrt(2.0*np.pi) * np.exp(-0.5*y**2)
    ax_res.plot(x, y, color='red')
    ax_res.minorticks_on()
    ax_res.xaxis.set_tick_params(labelbottom=True, labeltop=False)
    ax_res.set_xlim(0.0, 0.6)
    #===============================================

    pdf.savefig(fig)

  pdf.close()

  return

def plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend, flagnegresp, 
                 typetf, typemodel, resp_input, doshow=True, tf_lag_range=None, hist_lag_range=None, 
                 hist_bins=None, show_pmax=False, show_gap=False, labels=None):

  fp = open(fdir+fname)
  # read numbe of datasets
  line = fp.readline()
  nd = int(line[1:])

  # read number of data points in each dataset
  nl = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = np.array([int(i) for i in ls])
    for j in ns[1:]: # negelect the first one, which is number of continuum points
      if j == 0:
        typemodel = Model.nmap.value
    nl.append(ns)

  fp.close()

  if typemodel == Model.nmap.value: # nmap, only plot continuum
    plot_results_con(fdir, fname, flagtrend, doshow)
    return

  sample = np.atleast_2d(np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau))
  sample_info = np.loadtxt(fdir+"/data/posterior_sample_info1d.txt_%d"%ngau)
  idx_pmax = np.argmax(sample_info)
  data = np.loadtxt(fdir+fname)
  sall = np.loadtxt(fdir+"/data/pall.txt_%d"%ngau)

  if flagvar == 1:
    num_params_var = 3
  else:
    num_params_var = 3*nd
  
  if flagtrend > 0:
    trend = np.loadtxt(fdir+"/data/trend.txt_%d"%ngau)
  
  # number of parameters for long-term trend
  nq = flagtrend + 1

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
  
  # data sampling
  if show_gap is not None:
    if isinstance(show_gap, bool) and show_gap == True:
      DT_gap = []
      for i in range(nd):
        ns = nl[i]
        t_con = data[indx_con_data[i]:indx_con_data[i]+ns[0], 0]
        dt = t_con[1:] - t_con[:-1]
        span = t_con[-1] - t_con[0]
        ny = int(np.ceil(span/365.0))

        if ny > 1:
          dt = np.sort(dt)
          idx_dt = np.where(dt<365)[0]
          dt = dt[idx_dt]
          gap = np.mean(dt[-ny:])
        else:
          gap = None

        DT_gap.append(gap)

  # print time lags, median, and 68.3% confidence limits
  if flagnegresp == False:
    print("========No. of components: %d========="%ngau)
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
        if typemodel in [0, 2, 3]:  # general, vmap model
          for k in range(ngau):
            tt = int(typetf[k])
            indx = indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3
            if flagnegresp == 0:  # no negative response
              sample_lag[:] +=  _get_time_lag(sample, indx, tt) * np.exp(sample[:, indx+0])
              weight_lag[:] +=  np.exp(sample[:, indx+0])
            else:
              sample_lag[:] +=  sample[:, indx+1]
              weight_lag[:] +=  1.0

        elif typemodel == 1: # pmap model
          k = 0
          tt = int(typetf[k])
          indx = indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3
          sample_lag[:] +=  _get_time_lag(sample, indx, tt) * np.exp(sample[:, indx+0])
          weight_lag[:] +=  np.exp(sample[:, indx+0])
          indx0 = indx
          for k in range(1, ngau):
            tt = int(typetf[k])
            indx = indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3
            sample_lag[:] +=  _get_time_lag(sample, indx, tt) * np.exp(sample[:, indx+0] + sample[:, indx0+0])
            weight_lag[:] +=  np.exp(sample[:, indx+0] + sample[:, indx0+0])
    
        lag, err1, err2 = np.quantile(sample_lag/weight_lag, q=(0.5, (1.0-0.683)/2.0, 1.0-(1.0-0.683)/2.0))
        err1 = lag-err1
        err2 = err2 - lag
        print("Line %d: %.3f -%.3f +%.3f"%(j, lag, err1, err2))
  
  # load input respon function
  if resp_input != None:
    tran_input = np.loadtxt(resp_input)

  dtau = tau_upp - tau_low 
  ntau = 1000
  tran = np.zeros((sample.shape[0], ntau))
  
  shift = 0.0

  # open pdf file
  pdf = PdfPages(fdir+"/data/fig_%d.pdf"%ngau)
  print("Plotting to %s."%(fdir+"/data/fig_%d.pdf"%ngau))
  
  # set x-axis coordinates of sub figures
  if flagnegresp == False:
    figlc_center = 0.01
    figlc_centroid = 0.17
  else:
    figlc_center = 0.13
    # figlc_centroid = 0.22
  
  idx_q = 0 # index for long-term trend parameters
  for m in range(nd):
    ns = nl[m]
    ns_rec = nl_rec[m]
    fig = plt.figure(figsize=(12, 4+2*(len(ns)-2)))

    #====================================================================================================#
    # plot continuum
    con0 = data[indx_con_data[m]:indx_con_data[m]+ns[0], :]
    sall_con0 = sall[indx_con_rec[m]:(indx_con_rec[m]+ns_rec[0]), :] 
    
    axheight = 0.8/(len(ns))
    ax = fig.add_axes((0.49, 0.95-axheight, 0.35, axheight))
    
    ax.errorbar(con0[:, 0]-shift, con0[:, 1], yerr=con0[:, 2], ls='none', color='b', zorder=10, marker='o', markersize=1.5, elinewidth=0.5)
    
    ax.plot(sall_con0[:, 0]-shift, sall_con0[:, 1], color='k')
    ax.fill_between(sall_con0[:, 0]-shift, y1=sall_con0[:, 1]-sall_con0[:, 2], y2=sall_con0[:, 1]+sall_con0[:, 2], color='darkgrey')
    
    #===============================================
    # residuals
    if typemodel != 2: # do not plot residuals for vmap model
      ax_res = fig.add_axes((0.89, 0.95-axheight, 0.05, axheight))
      con0_intp = np.interp(con0[:, 0], sall_con0[:, 0], sall_con0[:, 1])
      ax_res.hist((con0[:, 1]-con0_intp)/con0[:, 2], bins=20, orientation='horizontal', density=True, color='C0')
      ax_res.yaxis.set_tick_params(labelleft=False, labelright=True)
      ax_res.yaxis.set_label_position("right")
      ax_res.set_ylabel("Res./Err.", labelpad=-5)
      ax_res.set_ylim(-4, 4)
      ylim_res = ax_res.get_ylim()
      ylim_res_max = max(abs(ylim_res[0]), abs(ylim_res[1]))
      ax_res.set_ylim(-ylim_res_max, ylim_res_max)
      y = np.linspace(-ylim_res_max, ylim_res_max, 100)
      x = 1.0/np.sqrt(2.0*np.pi) * np.exp(-0.5*y**2)
      ax_res.plot(x, y, color='red')
      ax_res.minorticks_on()
      ax_res.xaxis.set_tick_params(labelbottom=False)
      ax_res.set_xlim(0.0, 0.6)
    #===============================================

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
    # ax.set_ylabel('Flux')
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
    
    if labels is not None and isinstance(labels, list):
      if labels[m] is not None and labels[m][0] is not None:
        ax.text(0.05, 0.8, labels[m][0], transform=ax.transAxes)
    #====================================================================================================#
    # now do plotting lines
    
    # first determine range 
    if hist_lag_range is None:
      tau1, tau2, tau1_cent, tau2_cent = _get_hist_timelag_range(sample, indx_line, typetf, ngau, ns, m)
    else:
      tau1 = float(hist_lag_range[0])
      tau2 = float(hist_lag_range[1])
      tau1_cent, tau2_cent = tau1, tau2 
    
    if typemodel == 1:
      ratio1, ratio2 = _get_ratio_range(sample, indx_line, ngau, ns, m)
    
    if tf_lag_range is None:
      tau1_tf, tau2_tf = _get_tf_timelag_range(sample, indx_line, typetf, ngau, ns, m)

      tau1_tf = np.min((tau_low, tau1_tf))
      tau2_tf = np.max((tau_upp, tau2_tf))
    else:
      tau1_tf = float(tf_lag_range[0])
      tau2_tf = float(tf_lag_range[1])

    for j in range(1, len(ns)):
      hb = data[indx_con_data[m] + np.sum(ns[:j]):indx_con_data[m] + np.sum(ns[:j+1]), :] 
      sall_hb = sall[(indx_con_rec[m] + np.sum(ns_rec[:j])):(indx_con_rec[m] + np.sum(ns_rec[:j+1])), :]
      
      ##############################
      # first line light curve
      ##############################
      ax = fig.add_axes((0.49, 0.95-(j+1)*axheight, 0.35, axheight))
      
      ax.errorbar(hb[:, 0]-shift, hb[:, 1], yerr=hb[:, 2], ls='none', zorder=10, color='b', marker='o', markersize=1.5, elinewidth=0.5)
      
      ax.plot(sall_hb[:, 0]-shift, sall_hb[:, 1], color='k')
      ax.fill_between(sall_hb[:, 0]-shift, y1=sall_hb[:, 1]-sall_hb[:, 2], y2=sall_hb[:, 1]+sall_hb[:, 2], color='darkgrey')
      
      #===============================================
      # residuals
      ax_res = fig.add_axes((0.89, 0.95-(j+1)*axheight, 0.05, axheight))
      hb_intp = np.interp(hb[:, 0], sall_hb[:, 0], sall_hb[:, 1])
      ax_res.hist((hb[:, 1]-hb_intp)/hb[:, 2], bins=20, orientation='horizontal', density=True, color='C0')
      ax_res.yaxis.set_tick_params(labelleft=False, labelright=True)
      ax_res.yaxis.set_label_position("right")
      ax_res.set_ylabel("Res./Err.", labelpad=-5)
      ax_res.set_ylim(-4, 4)
      ylim_res = ax_res.get_ylim()
      ylim_res_max = max(abs(ylim_res[0]), abs(ylim_res[1]))
      ax_res.set_ylim(-ylim_res_max, ylim_res_max)
      y = np.linspace(-ylim_res_max, ylim_res_max, 100)
      x = 1.0/np.sqrt(2.0*np.pi) * np.exp(-0.5*y**2)
      ax_res.plot(x, y, color='red')
      ax_res.minorticks_on()
      ax_res.set_xlim(0.0, 0.6)
      if j != len(ns)-1:
        ax_res.xaxis.set_tick_params(labelbottom=False)
      else:
        ax_res.set_xlabel("Hist.")
      #===============================================

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
      # ax.set_ylabel("Flux")
      ax.set_xlim(xlim0[0], xlim0[1])
      
      # set ylim
      ymin = np.min(hb[:, 1])
      ymax = np.max(hb[:, 1])
      dy = ymax - ymin
      ax.set_ylim(ymin-0.1*dy, ymax+0.1*dy)

      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      ax.minorticks_on()
      if labels is not None and isinstance(labels, list):
        if labels[m] is not None and labels[m][j] is not None:
          ax.text(0.05, 0.8, labels[m][j], transform=ax.transAxes)
      ############################################
      # then posteriror distribution of time lags
      ############################################
      
      ax = fig.add_axes((figlc_center, 0.95-(j+1)*axheight, 0.15, axheight))

      for k in range(ngau):
        tt = int(typetf[k])
        if tt in [0, 1]: # gaussian or tophat
          cen = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          cen_pmax = sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        elif tt == 2:  # gamma, use peaks
          cen =  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          cen_pmax =  sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                +np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
        elif tt == 3:  # exp, use peaks
          cen =  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          cen_pmax =  sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        
        if hist_bins is None:
          cen_min = np.min(cen)
          cen_max = np.max(cen)
          bins = np.max((5, int((tau2-tau1)/(cen_max-cen_min + 1.0e-100) * 20)))
          bins = np.min((bins, 100))
        else:
          bins = hist_bins

        if k == 0:
          ax.hist(cen, density=True, range=(tau1, tau2), bins=bins, alpha=1)
        else:
          ax.hist(cen, density=True, range=(tau1, tau2), bins=bins, alpha=0.6)
        
        if show_pmax == True:
          ax.axvline(x=cen_pmax, ls='--', color='r')

      ax.set_xlim((tau1, tau2))
      
      # vmap model, no need to plot for the first lc, which has a zero lag wrt the driving lc. 
      if typemodel == 2 and j == 1 and ngau == 1:
        ax.set_visible(False)
      
      if flagnegresp == False:
        ax.yaxis.set_tick_params(labelleft=False)

      if (typemodel != 2 and j == 1) or (typemodel == 2 and j == 2 and ngau == 1 ) or (typemodel == 2 and j == 1 and ngau != 1):
        if typetf in [0, 1]:
          ax.set_title("Centers")
        else:
          ax.set_title("Peaks")

      ax.minorticks_on()
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")
    
      #====================================================================
      # only plot centroid lag when no negative response
      if flagnegresp == False:
        ax = fig.add_axes((figlc_centroid, 0.95-(j+1)*axheight, 0.15, axheight))
        
        if typemodel in [0, 2, 3]: # centroid time lag
          cent = np.zeros(sample.shape[0])
          norm = np.zeros(sample.shape[0])
          cent_pmax = 0.0
          norm_pmax = 0.0
          for k in range(ngau):
            tt = int(typetf[k])
            if tt in [0, 1]:
              norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                      * sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              
              norm_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                      * sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
              
            elif tt == 2:
              norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                    * (sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                      +2*np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
              
              norm_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                    * (sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                      +2*np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
            elif tt == 3:
              norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                    * (sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                      +np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
              
              norm_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
              cent_pmax += np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                    * (sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] \
                      +np.exp(sample[idx_pmax, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2]))
                        
          if hist_bins is None:
            cent_min = np.min(cent/norm)
            cent_max = np.max(cent/norm)
            bins = np.max((5, int((tau2_cent-tau1_cent)/(cent_max-cent_min + 1.0e-100) * 20)))
            bins = np.min((bins, 100))
          else:
            bins = hist_bins

          ax.hist(cent/norm, density=True, range=(tau1_cent, tau2_cent), bins=bins)
          ax.set_xlim((tau1_cent, tau2_cent))
          ax.minorticks_on()
          ax.yaxis.set_tick_params(labelleft=False)

          if show_pmax == True:
            ax.axvline(x = cent_pmax/norm_pmax, ls='--', color='r')

          if (typemodel != 2 and j == 1) or (typemodel == 2 and j == 2 and ngau == 1) or (typemodel==2 and j == 1 and ngau != 1):
            ax.set_title("Centroid")
            
          if j != len(ns)-1:
            ax.xaxis.set_tick_params(labelbottom=False)
          else:
            ax.set_xlabel("Time Lag (day)")

        elif typemodel == 1: # ratio hist
          
          for k in range(1, ngau):
            ratio = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3]/np.log(10.0)
            if k == 1:
              ax.hist(ratio, density=True, bins=30, alpha=1, range=(ratio1, ratio2))
            else:
              ax.hist(ratio, density=True, bins=30, alpha=0.6, range=(ratio1, ratio2))

            ax.yaxis.set_tick_params(labelleft=False)
            ax.minorticks_on()

            if j != len(ns)-1:
              ax.xaxis.set_tick_params(labelbottom=False)
            else:
              ax.set_xlabel(r"$\log R$")

            if j == 1:
              ax.set_title("Response Ratio")
        
        # vmap model, no need to plot for the first lc, which has a zero lag wrt the driving lc. 
        if typemodel == 2 and j == 1 and ngau == 1:
          ax.set_visible(False)
      
      ############################################
      # then transfer functions
      ############################################
      #===========================================================================================================
      # transfer function
      ax = fig.add_axes((0.33, 0.95-(j+1)*axheight, 0.15, axheight))
      
      tau = np.linspace(tau1_tf, tau2_tf, ntau)
      tran[:, :] = 0.0
      for i in range(sample.shape[0]):
        tran[i, :] = _calculate_tran(tau, sample[i, :], typemodel, typetf, ngau, flagnegresp, indx_line, m, j)
      
      tran_best = np.percentile(tran, 50.0, axis=0)
      tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
      tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)

      out = np.column_stack((tau, tran_best, tran_best-tran1, tran2-tran_best))
      np.savetxt(fdir+f"/data/tranfunc_{m}_{j}.txt_{ngau}", out, fmt="%f")
      
      ax.plot(tau, tran_best, color='k')
      ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')
      
      if show_pmax == True:
        tran_pmax = _calculate_tran(tau, sample[idx_pmax, :], typemodel, typetf, ngau, flagnegresp, indx_line, m, j)
        ax.plot(tau, tran_pmax, label=r'$L_{\rm max}$', color='r', ls='--')
        

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
      
      if resp_input != None or show_pmax == True:
        ax.legend(fontsize=10)

      ax.set_xlim(tau1_tf, tau2_tf)

      ylim = ax.get_ylim()
      if resp_input == None:
        ymax = np.max(tran_best)
        ymin = np.min(tran_best)
        dy = ymax - ymin 
        ax.set_ylim(np.max((ylim[0], ymin-0.1*dy)), np.min((ylim[1], np.max(np.max(tran_best)*1.5))))
      else:
        ymax = np.max((np.max(tran_best), np.max(tran_input[:, 1])))
        ymin = np.min((np.min(tran_best), np.min(tran_input[:, 1])))
        dy = ymax - ymin 
        ax.set_ylim(np.max((ylim[0], ymin-0.1*dy)), np.min((ylim[1], np.max((np.max(tran_best)*1.5, np.max(tran_input[:, 1])*1.5)))))
      
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
      
      if show_gap is not None:
        if  isinstance(show_gap, bool) and show_gap == True:
          xlim = ax.get_xlim()
          gap = DT_gap[m]
          offset = 0
          if gap is not None:
            while gap + offset > xlim[0] and gap + offset < xlim[1]:
              ylim = ax.get_ylim()
              ax.fill_between(x=[offset+365/2-gap/2, offset+365/2+gap/2], y1=[ylim[1], ylim[1]], y2=[ylim[0], ylim[0]], color='darkgrey', alpha=0.5)
              ax.set_ylim(ylim[0], ylim[1])
              ax.text(offset+365/2, ylim[1]-0.1*(ylim[1]-ylim[0]), "gap", ha='center', fontsize=10)
              offset += 365
          ax.set_xlim(xlim[0], xlim[1])

        elif type(show_gap) == list or type(show_gap) == np.ndarray:
          xlim = ax.get_xlim()
          center = float(show_gap[m*2])
          gap = float(show_gap[m*2+1])
          offset = 0
          while gap + offset > xlim[0] and gap + offset < xlim[1]:
            ylim = ax.get_ylim()
            ax.fill_between(x=[offset+center-gap/2, offset+center+gap/2], y1=[ylim[1], ylim[1]], y2=[ylim[0], ylim[0]], color='darkgrey', alpha=0.5)
            ax.set_ylim(ylim[0], ylim[1])
            ax.text(offset+center, ylim[1]-0.1*(ylim[1]-ylim[0]), "gap", ha='center', fontsize=10)
            offset += 365
          ax.set_xlim(xlim[0], xlim[1])
      
    pdf.savefig(fig)

    if doshow:
      plt.show()
    else:
      plt.close()

    idx_q += len(ns) * nq

  pdf.close()
  return

def plot_results_all(args, param, doshow=True, tf_lag_range=None, hist_lag_range=None, hist_bins=None):
  """
  plot function called by directly execuating "plotfig.py" in terminal. 
  """
  try:
    fdir = param["FileDir"]+"/"
  except:
    raise IOError("FileDir is not set!")
  
  try:
    fname = param["DataFile"]
  except:
    raise IOError("DataFile is not set!")
  
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
  
  typetf = None
  typemodel = int(param["TypeModel"])
  if typemodel != 3:
    try:  
      ngau_low = int(param["NumCompLow"])
    except:
      raise IOError("NumCompLow is not set!")

    try: 
      ngau_upp = int(param["NumCompUpp"])
    except:
      raise IOError("NumCompUpp is not set!")
  else: # mmap mode
    typetf = param["StrTypeTFMix"]
    ngau_low = len(typetf)
    ngau_upp = ngau_low
  
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
    flagnegresp = int(param["FlagNegativeResp"])
  except:
    flagnegresp = 0
  
  show_gap = False
  if args.show_gap:
    try:
      str_gap = param["StrGapPrior"]
      show_gap = str_gap[1:-1].split(":")
    except:
      show_gap = args.show_gap
  
  for ngau in range(ngau_low, ngau_upp+1):
    if typemodel != 3:
      typetf = param["TypeTF"]*ngau
    
    plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend, 
                 flagnegresp, typetf, typemodel, args.resp_input, 
                 doshow=doshow, tf_lag_range=args.tf_lag_range, 
                 hist_lag_range=args.hist_lag_range, hist_bins=args.hist_bins, 
                 show_pmax=args.show_pmax, show_gap=show_gap)

def _param_parser(fname):
  """
  parse parameter file
  """
  config = cp.RawConfigParser(delimiters=' ', comment_prefixes='#', inline_comment_prefixes='#', 
  default_section=cp.DEFAULTSECT, empty_lines_in_values=False, allow_no_value=True)
  
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
  parser.add_argument('--hist_bins', type=int, nargs='+', help="number of bins for the histograms, e.g., --hist_bins 20")
  parser.add_argument('--show_gap', action='store_true', default=False, help="whether show seasonal gaps, e.g., --show_gap")
  parser.add_argument('--show_pmax', action='store_true', default=False, help="whether show the results of the maximum posterior ppint, e.g., --show_pmax")
  args = parser.parse_args()

  if args.param == None:
    print("Please specify paramter file!")
    print("e.g., python plotfig.py --param src/param")
    print(parser.parse_args(['-h']))
    sys.exit()

  fparam = args.param
  param = _param_parser(fparam)
  
  if "TypeModel" not in param:
    param["TypeModel"] = "0"

  plot_results_all(args, param)
# MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
#

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import (MultipleLocator)
import sys, os
import configparser as cp 
from matplotlib.backends.backend_pdf import PdfPages
import argparse

# set the default parameters for matplotlib
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot_line_decomp(fdir, fname, ngau, tau_low, tau_upp, typetf, typemodel, flagnegresp, resp_input, doshow=True):

  data = np.loadtxt(fdir+fname)
  pall = np.loadtxt(fdir+"data/pall.txt_%d"%ngau)
  
  #red number of data points
  fp = open(fdir+fname, "r")
  line = fp.readline()
  nd = int(line[1:])

  nls_data = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = np.array([int(i) for i in ls])
    nls_data.append(ns)
  fp.close()
  
  #red number of points
  fp = open(fdir+"data/pall.txt_%d"%ngau, "r")
  line = fp.readline()

  nls = []
  for i in range(nd):
    line = fp.readline()
    ls = line[1:].split(":")
    ns = np.array([int(i) for i in ls])
    nls.append(ns)
  fp.close()
  
  # load component
  comps = []
  for i in range(ngau):
    comps.append(np.loadtxt(fdir+"data/pline.txt_%d_comp%d"%(ngau, i)))
  
  # load posterior samples
  sample = np.loadtxt(fdir+"data/posterior_sample1d.txt_%d"%ngau)
  
  # load input respon function
  if resp_input != None:
    tran_input = np.loadtxt(resp_input)
  
  if flagnegresp==True:
    fq = np.loadtxt(fdir+"data/trend.txt_%d"%ngau)

  # determine the maximum and minimum delays
  tau_min =  1.0e10
  tau_max = -1.0e10
  idx = nd*3  # variability parameters come first, each dataset has 3 parameters
  for i in range(nd):
    didx = (len(nls_data[i])-1)*(ngau*3+1) # each line has (ngau*3+1) parameters
    tau_min = min(tau_min, np.min(sample[:, idx+1+1:idx+didx:3]-np.exp(sample[:, idx+1+2:idx+didx:3])))
    tau_max = max(tau_max, np.max(sample[:, idx+1+1:idx+didx:3]+np.exp(sample[:, idx+1+2:idx+didx:3])))

    idx += didx
  
  tau_min = np.min((tau_low, tau_min))
  tau_max = np.max((tau_max, tau_upp))

  if resp_input is not None:
    tau_min = np.min((tau_min, tran_input[0, 0]))
    tau_max = np.max((tau_max, tran_input[-1, 0]))

  tau = np.linspace(tau_min, tau_max, 1000)
  tran = np.zeros((sample.shape[0], 1000))
  
  pdf = PdfPages(fdir+"data/fig_line_decomp_%d.pdf"%ngau)
  
  nq = 1
  idx_q = 0 # index for long-term trend parameters
  for m in range(nd):
    nl_data = nls_data[m]
    nl = nls[m]
    con_data = data[:nl_data[0], :]
    con_rec = pall[:nl[0], :]
    
    idx_hb_data = nl_data[0]
    idx_hb = nl[0]

    tran[:, :] = 0.0
    for j in range(1, len(nl_data)):

      hb_data = data[idx_hb_data:idx_hb_data+nl_data[j], :]
      hb_rec = pall[idx_hb:idx_hb+nl[j], :]
      if typetf == 0:
        for i in range(sample.shape[0]):
          # loop over gaussians
          if typemodel == 0:  # general model
            for k in range(ngau):

              if flagnegresp == False:
                amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)

          elif typemodel == 1:  # pmap model
            k = 0
            amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
            tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)
            for k in range(1, ngau):
              amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)
      elif typetf == 1:
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0:   # general model
            for k in range(ngau):

              if flagnegresp == False:
                amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))

          elif typemodel == 1: # pmap model
            k = 0
            amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
            tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))
            for k in range(1, ngau):
              amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, :] += amp/sig/2.0 *(np.heaviside(sig-np.abs(tau-cen), 1.0))
      elif typetf == 2:
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0 or typemodel == 2: # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              idx_tau = np.where(tau >= cen)[0]
              tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)

          elif typemodel == 1: #pmap model 
            k = 0
            amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
            idx_tau = np.where(tau >= cen)[0]
            tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
            for k in range(1, ngau):
              amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
      else:  # gamma
        for i in range(sample.shape[0]):
          # loop over tophats
          if typemodel == 0 or typemodel == 2: # general, vmap model
            for k in range(ngau):

              if flagnegresp == 0:
                amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
              else:
                amp =      sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0]

              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              
              idx_tau = np.where(tau >= cen)[0]
              tran[i, idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)

          elif typemodel == 1: #pmap model 
            k = 0
            amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0])
            cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
            sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
            idx_tau = np.where(tau >= cen)[0]
            tran[i, idx_tau] += amp/sig**2 * (tau[idx_tau]-cen) * np.exp(-(tau[idx_tau]-cen)/sig)
            for k in range(1, ngau):
              amp = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+0] + \
                           sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+0*3+0])
              
              cen =        sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+1]
              sig = np.exp(sample[i, 3*nd + (j-1)*(ngau*3+1) + 1+k*3+2])
              tran[i, idx_tau] += amp/sig * np.exp(-(tau[idx_tau]-cen)/sig)

     
      tran_best = np.percentile(tran, 50.0, axis=0)
      tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
      tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)

      fig = plt.figure(figsize=(12, 9))
      ax = fig.add_axes((0.1, 0.8, 0.8, 0.19))
      ax.plot(tau, tran_best, color='k')
      ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')
      ax.minorticks_on()

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

      ax.set_xlabel("Time Lag")
      ax.set_ylabel("Transfer Function")
      ylim = ax.get_ylim()
      ax.set_ylim(ylim[0], np.min((ylim[1], np.max(np.max(tran_best)*1.5))))
      
      # plot light curves
      ax = fig.add_axes((0.1, 0.51, 0.8, 0.19))
      ax.errorbar(con_data[:, 0], con_data[:, 1], yerr=con_data[:, 2], ls='none', marker='o', markersize=3, color='k', markerfacecolor='C0', markeredgewidth=0.4, elinewidth=0.8)
      ax.plot(pall[:nl[0], 0], pall[:nl[0], 1], lw=1)
      ax.fill_between(pall[:nl[0],0], y1=pall[:nl[0], 1]-pall[:nl[0], 2], y2=pall[:nl[0], 1]+pall[:nl[0], 2], color='darkgrey', zorder=0)
      ax.set_ylabel(r"$F_\lambda$(cont)")
      #ax.set_xlim((-10.0, 260.0))
      ax.set_xticklabels([])
      ax.minorticks_on()
      
      line_error = 0
      ax = fig.add_axes((0.1, 0.1, 0.8, 0.40))
      ax.errorbar(hb_data[:, 0], hb_data[:, 1], yerr=np.sqrt(hb_data[:, 2]**2+line_error**2), ls='none', marker='o', markersize=3, color='k', markerfacecolor='C0', markeredgewidth=0.4, elinewidth=0.8)
      
      # plot mean of line
      if flagnegresp == True:
        mean = fq[idx_q+1, 0]
        ax.axhline(y=mean, ls='--', lw=1, color='grey')
      
      for i in range(ngau):
        l0 = comps[i][idx_hb:idx_hb+nl[j], :]
        lg0,=ax.plot(l0[:, 0], l0[:, 1], lw=1, label='%d'%(i+1), color="C%d"%i)
        ax.fill_between(l0[:, 0], y1=l0[:, 1]-l0[:, 2], y2=l0[:, 1]+l0[:, 2], color="C%d"%i, alpha=0.5, zorder=0)

      lgt,=ax.plot(hb_rec[:, 0], hb_rec[:, 1], lw=1, label='total')
      ax.fill_between(hb_rec[:, 0], y1=hb_rec[:, 1]-hb_rec[:, 2], y2=hb_rec[:, 1]+hb_rec[:, 2], color='darkgrey', zorder=0)

      ax.set_xlabel(r"Time")
      ax.set_ylabel(r"$F(\rm line$)")
      #ax.set_xlim((-10.0, 260.0))
      ylim = ax.get_ylim()
      ax.legend()
      ax.minorticks_on()

      if doshow==True:
        plt.show()
      else:
        plt.close()

      pdf.savefig(fig)

      idx_hb_data += nl_data[j]
      idx_hb += nl[j]
    
    idx_q = len(nl_data) * nq

  pdf.close()

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
  parser = argparse.ArgumentParser(usage="python plotdecomp.py [options]")
  parser.add_argument('--param', type=str, help="parameter file")
  parser.add_argument('--resp_input', type=str, help="str, a file storing input response function")
  parser.add_argument('--tf_lag_range', type=float, nargs='+', help="time lag range for the transfer function, e.g., --tf_lag_range 0 100")
  parser.add_argument('--hist_lag_range', type=float, nargs='+', help="time lag range for the histograms, e.g., --hist_lag_range 0 100")
  args = parser.parse_args()

  if args.param == None:
    print("Please specify paramter file!")
    print("e.g., python plotdecomp.py --param src/param")
    print(parser.parse_args(['-h']))
    sys.exit(0)

  fparam = args.param
  param = _param_parser(fparam)

  # default option
  if "TypeModel" not in param:
    param["TypeModel"] = "0"
  
  try:
    fdir = param["FileDir"]+"/"
  except:
    raise IOError("FileDir is not set!")

  try:
    fname = param["DataFile"]
  except:
    raise IOError("DataFile is not set!")

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
  
  if ngau_upp < 2:
    print("The number of Gaussian < 2, no need to show decomposition!")
    sys.exit()

  for ngau in range(2, ngau_upp+1):
    plot_line_decomp(fdir, fname, ngau, tau_low, tau_upp, typetf, typemodel, flagnegresp, args.resp_input)


import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os
import configparser as cp 


def plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend):
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif', size=15)

  sample = np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau)
  data = np.loadtxt(fdir+fname)
  sall = np.loadtxt(fdir+"/data/pall.txt_%d"%ngau)
  scale = int(sall.shape[0]/data.shape[0])

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
  
  # assign index of cont data
  indx_con_data = []
  indx_con_data.append(0)
  for i in range(1, nd):
    ns = nl[i-1]
    indx_con_data.append(np.sum(ns) + indx_con_data[i-1])

  # assign index of the parmaeter for the first line of each dataset 
  indx_line = []
  indx_line.append(num_params_var)
  for i in range(1, nd):
    if flagtran == 1:
      indx_line.append(num_params_var)
    else:
      indx_line.append(indx_line[i-1] + (len(nl[i-1])-1)*(1+ngau*3))

  # print time lags, median, and 68.3% confidence limits
  print("========No. Gaussian: %d========="%ngau)
  print("ID:      lag     -elo      +eup")
  sample_lag = np.zeros(sample.shape[0])
  weight_lag = np.zeros(sample.shape[0])
  #loop over datasets
  for m in range(nd):
    print("Dataset %d"%m)
    #loop over lines
    for j in range(1, len(ns)):
      sample_lag[:] = 0.0
      weight_lag[:] = 0.0
      for k in range(ngau):
        sample_lag[:] +=  sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1] * np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        weight_lag[:] +=  np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
  
      lag, err1, err2 = np.quantile(sample_lag/weight_lag, q=(0.5, (1.0-0.683)/2.0, 1.0-(1.0-0.683)/2.0))
      err1 = lag-err1
      err2 = err2 - lag
      print("Line %d: %.3f -%.3f +%.3f"%(j, lag, err1, err2))

  dtau = tau_upp - tau_low 
  ntau = 500
  tran = np.zeros((sample.shape[0], ntau))
  
  shift = 0.0
  
  # open pdf file
  pdf = PdfPages(fdir+"/data/fig_%d.pdf"%ngau)

  idx_q = 0 # index for long-term trend parameters
  for m in range(nd):
    ns = nl[m]
    fig = plt.figure(figsize=(12, 4+2*(len(ns)-2)))
    
    #===================================
    # plot continuum
    con0 = data[indx_con_data[m]:indx_con_data[m]+ns[0], :]
    sall_con0 = sall[indx_con_data[m]*scale:(indx_con_data[m]+ns[0])*scale, :] 
    
    axheight = 0.8/(len(ns))
    ax = fig.add_axes((0.56, 0.95-axheight, 0.35, axheight))
    
    ax.errorbar(con0[:, 0]-shift, con0[:, 1], yerr=con0[:, 2], ls='none', color='b', zorder=10)
    
    ax.plot(sall_con0[:, 0]-shift, sall_con0[:, 1], color='k', lw=1)
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
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xlim0 = xlim 
  
    # plot line
    # set time lag range for Gaussian centers and centriods
    tau1 = 1.0e10
    tau2 = -1.0e10
    for j in range(1, len(ns)):      
      for k in range(ngau):
        tau1 = np.min((tau1, np.min(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))
        tau2 = np.max((tau2, np.max(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1])))

    # set time lag range for transfer function 
    tau1_tf = 1.0e10
    tau2_tf = -1.0e10
    for j in range(1, len(ns)):   
      for k in range(ngau):
        tau1_tf = np.min((tau1_tf, np.min(sample[:, indx_line[m]+1+k*3+1]-2*np.exp(sample[:, indx_line[m]+1+k*3+2]))))
        tau2_tf = np.max((tau2_tf, np.max(sample[:, indx_line[m]+1+k*3+1]+2*np.exp(sample[:, indx_line[m]+1+k*3+2]))))
      
    tau1_tf = np.min((tau_low, tau1_tf))
    tau2_tf = np.max((tau_upp, tau2_tf))

    # now do plotting
    for j in range(1, len(ns)):
      hb = data[indx_con_data[m] + np.sum(ns[:j]):indx_con_data[m] + np.sum(ns[:j+1]), :] 
      sall_hb = sall[(indx_con_data[m] + np.sum(ns[:j]))*scale:(indx_con_data[m] + np.sum(ns[:j+1]))*scale, :]
      
      # histogram of time lags 
      ax = fig.add_axes((0.05, 0.95-(j+1)*axheight, 0.16, axheight))
      for k in range(ngau):
        cen = sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
        if k == 0:
          ax.hist(cen, density=True, range=(tau1, tau2), bins=20, alpha=1)
        else:
          ax.hist(cen, density=True, range=(tau1, tau2), bins=20, alpha=0.6)

      ax.set_xlim((tau1, tau2))
      ax.yaxis.set_tick_params(labelleft=False)
      if j == 1:
        ax.set_title("Gaussian Centers")
      ax.minorticks_on()
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")

      # centroid time lag
      ax = fig.add_axes((0.22, 0.95-(j+1)*axheight, 0.16, axheight))

      cent = np.zeros(sample.shape[0])
      norm = np.zeros(sample.shape[0])
      for k in range(ngau):
        norm += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
        cent += np.exp(sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0]) \
                * sample[:, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
      
      ax.hist(cent/norm, density=True, range=(tau1, tau2), bins=20)
      ax.set_xlim((tau1, tau2))
      ax.minorticks_on()
      ax.yaxis.set_tick_params(labelleft=False)
      if j == 1:
        ax.set_title("Centroid Time Lag")
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")

      # transfer function
      ax = fig.add_axes((0.39, 0.95-(j+1)*axheight, 0.16, axheight))
      
      tau = np.linspace(tau1_tf, tau2_tf, ntau)
      tran[:, :] = 0.0
      for i in range(sample.shape[0]):
        # loop over gaussians
        for k in range(ngau):
          amp = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+0])
          cen =        sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+1]
          sig = np.exp(sample[i, indx_line[m] + (j-1)*(ngau*3+1) + 1+k*3+2])
          tran[i, :] += amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)
      
      tran_best = np.percentile(tran, 50.0, axis=0)
      tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
      tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)
      
      ax.plot(tau, tran_best, color='k')
      ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')
      ax.set_xlim((tau[0], tau[-1]))
      ax.set_ylim(0.0, np.min((ylim[1], np.max(tran_best)*1.5)))
      
      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("Time Lag (day)")
      
      if j == 1:
        ax.set_title("Transfer Function")
    
      ax.yaxis.set_tick_params(labelleft=False)
      ax.minorticks_on()
      #if(tau[0]<0.0):
      #  ax.axvline(x=0.0, ls='--', color='red')
      
      # then line light curve
      ax = fig.add_axes((0.56, 0.95-(j+1)*axheight, 0.35, axheight))
      
      ax.errorbar(hb[:, 0]-shift, hb[:, 1], yerr=hb[:, 2], ls='none', zorder=10, color='b')
      
      ax.plot(sall_hb[:, 0]-shift, sall_hb[:, 1], color='k')
      ax.fill_between(sall_hb[:, 0]-shift, y1=sall_hb[:, 1]-sall_hb[:, 2], y2=sall_hb[:, 1]+sall_hb[:, 2], color='darkgrey')
      
      # plot long-term trend 
      if flagtrend > 0:
        xlim = ax.get_xlim()
        x = np.linspace(sall_hb[0, 0], sall_hb[-1, 0], 100)
        y = np.zeros(100)
        for k in range(nq):
          y+= trend[idx_q + 1*nq + k, 0]* x**(k)
        
        ax.plot(x, y, ls='--', color='grey')

      if j != len(ns)-1:
        ax.xaxis.set_tick_params(labelbottom=False)
      else:
        ax.set_xlabel("JD")
      ax.yaxis.set_tick_params(labelleft=False)
      ax.yaxis.set_tick_params(labelright=True)
      ax.yaxis.set_label_position("right")
      ax.set_ylabel("Flux")
      ax.set_xlim(xlim[0], xlim[1])
      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      ax.minorticks_on()

      idx_q += len(ns) * nq
  
    plt.show()
    pdf.savefig(fig)
  
  pdf.close()
  return

def plot_results_all(fdir, fname, ngau_low, ngau_upp, tau_low, tau_upp, flagvar, flagtran, flagtrend):
  for ngau in range(ngau_low, ngau_upp+1):
    plot_results(fdir, fname, ngau, tau_low, tau_upp, flagvar, flagtran, flagtrend)

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
  if(len(sys.argv) < 2):
    print("Please specify paramter file!")
    print("e.g., python plotfig.py src/param")
    exit(0)
  fparam = sys.argv[1]
  param = _param_parser(fparam)
  fdir = param["FileDir"]
  flagvar = int(param["FlagUniformVarParams"])
  flagtran = int(param["FlagUniformTranFuns"])
  flagtrend = int(param["FlagLongtermTrend"])
  ngau_low = int(param["NumGaussianLow"])
  ngau_upp = int(param["NumGaussianUpp"])
  tau_low = float(param["LagLimitLow"])
  tau_upp = float(param["LagLimitUpp"])
  fname = param["DataFile"]
  plot_results_all(fdir, fname, ngau_low, ngau_upp, tau_low, tau_upp, flagvar, flagtran, flagtrend)

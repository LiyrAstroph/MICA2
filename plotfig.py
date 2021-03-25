import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os
import configparser as cp 


def plot_results(fdir, fname, ngau, tau_low, tau_upp):
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif', size=18)

  sample = np.loadtxt(fdir+"/data/posterior_sample1d.txt_%d"%ngau)
  data = np.loadtxt(fdir+fname)
  sall = np.loadtxt(fdir+"/data/pall.txt_%d"%ngau)
  scale = int(sall.shape[0]/data.shape[0])

  fp = open(fdir+fname)
  line = fp.readline()
  line = fp.readline()
  ls = line[1:].split(":")
  ns = np.array([int(i) for i in ls])
  fp.close()
  
  dtau = tau_upp - tau_low 
  tau = np.linspace(tau_low-0.5*dtau, tau_upp+0.5*dtau, 100)
  tran = np.zeros((sample.shape[0], 100))
  
  shift = 0.0

  fig = plt.figure(figsize=(8, 4+2*(len(ns)-2)))

  #===================================
  con0 = data[:np.sum(ns[:1]), :]
  sall_con0 = sall[:np.sum(ns[:1])*scale, :] 
  
  axheight = 0.8/(len(ns))
  ax = fig.add_axes((0.4, 0.95-axheight, 0.5, axheight))
  
  ax.errorbar(con0[:, 0]-shift, con0[:, 1], yerr=con0[:, 2], ls='none', color='b', zorder=10)
  
  ax.plot(sall_con0[:, 0]-shift, sall_con0[:, 1], color='k', lw=1)
  ax.fill_between(sall_con0[:, 0]-shift, y1=sall_con0[:, 1]-sall_con0[:, 2], y2=sall_con0[:, 1]+sall_con0[:, 2], color='darkgrey')
  
  ax.xaxis.set_tick_params(labeltop=False)
  ax.xaxis.set_tick_params(labelbottom=False)
  ax.yaxis.set_tick_params(labelleft=False)
  ax.yaxis.set_tick_params(labelright=True)
  ax.yaxis.set_label_position("right")
  ax.set_ylabel('Flux')
  #ax.xaxis.set_label_position("top")
  #ax.set_xlabel('HJD - 2450000')
  
  xlim = ax.get_xlim()
  ylim = ax.get_ylim()
  xlim0 = xlim 


  for j in range(1, len(ns)):
    hb = data[np.sum(ns[:j]):np.sum(ns[:j+1]), :] 
    sall_hb = sall[np.sum(ns[:j])*scale:np.sum(ns[:j+1])*scale, :] 
    ax = fig.add_axes((0.08, 0.95-(j+1)*axheight, 0.3, axheight))
    
    for i in range(sample.shape[0]):
      amp = np.exp(sample[i, 3+(j-1)*4+1])
      cen =        sample[i, 3+(j-1)*4+2]
      sig = np.exp(sample[i, 3+(j-1)*4+3])
      tran[i, :] = amp/sig * np.exp(-0.5*(tau - cen)**2/sig**2)
    
    tran_best = np.percentile(tran, 50.0, axis=0)
    tran1 = np.percentile(tran, (100.0-68.3)/2.0, axis=0)
    tran2 = np.percentile(tran, 100.0-(100.0-68.3)/2.0, axis=0)
    
    ax.plot(tau, tran_best, color='k')
    ax.fill_between(tau, y1=tran1, y2=tran2, color='darkgrey')
    ax.set_xlim((tau[0], tau[-1]))
    
    if j != len(ns)-1:
      ax.xaxis.set_tick_params(labelbottom=False)
    else:
      ax.set_xlabel("Time Lag (day)")
  
    ax.set_ylabel("TF")
  
    ax.yaxis.set_tick_params(labelleft=False)
    ax.minorticks_on()
    if(tau[0]<0.0):
      ax.axvline(x=0.0, ls='--', color='red')
    
    ax = fig.add_axes((0.4, 0.95-(j+1)*axheight, 0.5, axheight))
    
    ax.errorbar(hb[:, 0]-shift, hb[:, 1], yerr=hb[:, 2], ls='none', zorder=10, color='b')
    
    ax.plot(sall_hb[:, 0]-shift, sall_hb[:, 1], color='k')
    ax.fill_between(sall_hb[:, 0]-shift, y1=sall_hb[:, 1]-sall_hb[:, 2], y2=sall_hb[:, 1]+sall_hb[:, 2], color='darkgrey')
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

  plt.show()
  fig.savefig(fdir+"/data/fig_%d.pdf"%ngau, bbox_inches='tight')


def plot_results_all(fdir, fname, ngau_low, ngau_upp, tau_low, tau_upp):
  for ngau in range(ngau_low, ngau_upp+1):
    plot_results(fdir, fname, ngau, tau_low, tau_upp)

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
    exit(0)
  fparam = sys.argv[1]
  param = _param_parser(fparam)
  fdir = param["FileDir"]
  ngau_low = int(param["NumGaussianLow"])
  ngau_upp = int(param["NumGaussianUpp"])
  tau_low = float(param["LagLimitLow"])
  tau_upp = float(param["LagLimitUpp"])
  fname = param["DataFile"]
  plot_results_all(fdir, fname, ngau_low, ngau_upp, tau_low, tau_upp)

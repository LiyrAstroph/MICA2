# MICA2
# A code for time lag measurement in reverberation mapping
# 
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Jun 22, 2018
#

import numpy as np
import matplotlib.pyplot as plt 
from numpy import fft 
import copy
import sys

__all__ = ["simlc", "format_mica"]

#=======================================================
# DRW PSD function
#
#=======================================================
def psd_drw(fk, arg, freq_limit):
  sigma, tau, cnoise = arg[:3]
  
  nud = 1.0/(2.0*np.pi*tau)
  psd = np.zeros(fk.shape)
  
  psd[:] = 2.0*sigma**2*tau/(1.0 + (fk/nud)**2) + cnoise
    
  return psd

#=======================================================
# square root of DRW PSD function
#
#=======================================================
def psd_drw_sqrt(fk, arg, freq_limit):
  sigma, tau, cnoise = arg[:3]
  nud = 1.0/(2.0*np.pi*tau)
  psd = np.zeros(fk.shape)
  
  psd[:] = 2.0*sigma**2*tau/(1.0 + (fk/nud)**2) + cnoise
  
  return np.sqrt(psd)


#=======================================================
# single power law PSD function
#
#=======================================================
def psd_power_law(fk, arg, freq_limit):
  A, alpha, cnoise =arg[:3]
  
  psd = np.zeros(fk.shape)
  
  idx = np.where(fk >= freq_limit)
  psd[idx[0]] = A * fk[idx[0]]**(-alpha) + cnoise
  idx = np.where(fk < freq_limit)
  psd[idx[0]] = A * freq_limit**(-alpha) + cnoise
  
  return psd

#=======================================================
# square root of single power law PSD function
#
#=======================================================
def psd_power_law_sqrt(fk, arg, freq_limit):
  A, alpha, cnoise =arg[:3]
  
  psd = np.zeros(fk.shape)
  
  idx = np.where(fk >= freq_limit)
  psd[idx[0]] = A * fk[idx[0]]**(-alpha) + cnoise
  idx = np.where(fk < freq_limit)
  psd[idx[0]] = A * freq_limit**(-alpha) + cnoise
  
  return np.sqrt(psd)

#=======================================================  
# generate time series with given power-law PSD
#
# note that integrating PSD over (-infty, +infty) equals 
# to variation of light curves
#=======================================================
def genlc_psd_pow(model, nd, DT, freq_limit):  
  W = 4
  V = 4
  nd_sim = int(W*V*nd)
  arg = model
  fft_work = np.zeros(nd_sim//2+1, dtype=complex)
  
  fft_work[0] = np.random.randn()+1j*0.0
  
  freq = 1.0/(nd_sim * DT/W) * np.linspace(0.0, nd_sim//2, nd_sim//2+1)
  fft_work[1:nd_sim//2] = psd_power_law_sqrt(freq[1:nd_sim//2], arg, freq_limit)/np.sqrt(2.0) \
                        * ( np.random.randn(int(nd_sim//2-1)) + 1j*np.random.randn(int(nd_sim//2-1)) )
  fft_work[nd_sim//2:] = psd_power_law_sqrt(freq[nd_sim//2:], arg, freq_limit) * (np.random.randn() + 1j*0.0)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT/W * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0)
  fs = fs*norm
  
  return ts[::W], fs[::W]
  
#=======================================================  
# generate time series with given drw PSD
# 
# note that integrating PSD over (-infty, +infty) equals 
# to variation of light curves
#=======================================================
def genlc_psd_drw(model, nd, DT, freq_limit):  
  W = 10
  V = 10
  nd_sim = int(W*V*nd)
  arg = model
  fft_work = np.zeros(nd_sim//2+1, dtype=complex)
  
  fft_work[0] = np.random.randn()+1j*0.0
  
  freq = 1.0/(nd_sim * DT/W) * np.linspace(0.0, nd_sim//2, nd_sim//2+1)
  fft_work[1:nd_sim//2] = 2.0**0.5*psd_drw_sqrt(freq[1:nd_sim//2], arg, freq_limit)/np.sqrt(2.0) \
                        * ( np.random.randn(int(nd_sim//2-1)) + 1j*np.random.randn(int(nd_sim//2-1)) )
  fft_work[nd_sim//2:] = 2.0**0.5*psd_drw_sqrt(freq[nd_sim//2:], arg, freq_limit) * (np.random.randn() + 1j*0.0)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT/W * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0)
  fs = fs*norm
  
  return ts[::W], fs[::W]  

###################################################
# DRW
def drw(tot, model):  
  freq_limit = 1.0e-4
  DT = 0.01
  ts, fs = genlc_psd_drw(model, int(tot/DT), DT, freq_limit)
  fs += 1.0
  
  return ts, fs

##################################################
def psd(tot, model):
  freq_limit = 5.0e-3
  DT = 0.01
  ts, fs = genlc_psd_pow(model, int(tot/DT), DT, freq_limit)
  fs += 1.0
  
  return ts, fs
##################################################
def convolve_fft(con, resp):
  """
  convolution using FFT
  """
  resp_pad = np.zeros(con.shape[0])
  resp_pad[:resp.shape[0]] = resp
  con_fft = fft.rfft(con)
  resp_fft = fft.rfft(resp_pad)
  conv_fft = con_fft * resp_fft 
  conv = fft.irfft(conv_fft, n = con.shape[0])
  return conv
##################################################
def format_mica(fname, data1, data2):
  """
  This function generates a formatted file for MICA.
  If data1 is none, the vmap mode in MICA is presumed.

  Parameters
  ----------
  fname : string
    File name
  
  data1 : 2D array like
    The driving light curve, time, flux, and error.
  
  data2 : 2D array like
    The responding light curve, time, flux, and error.
  
  Returns
  -------
  None : None
    No returns.
  """
  
  # data1 is none
  if data1 is None:
    if not isinstance(data2, list):
      fp = open(fname, "w")
      fp.write("# 1\n")
      fp.write("# %d:%d\n"%(0, data2.shape[0]))
      np.savetxt(fp, data2, fmt="%f")
      fp.close()
    else:
      fp = open(fname, "w")
      fp.write("# 1\n")
      fp.write("# %d"%0)
      for i in range(len(data2)):
        fp.write(":%d"%data2[i].shape[0])
      fp.write("\n")
      for i in range(len(data2)):
        np.savetxt(fp, data2[i], fmt="%f")
        fp.write("\n")
      
      fp.close()

  else:
    # data1 has only one set
    if not isinstance(data1, list):
      if not isinstance(data2, list):
        fp = open(fname, "w")
        fp.write("# 1\n")
        fp.write("# %d:%d\n"%(data1.shape[0], data2.shape[0]))
        np.savetxt(fp, data1, fmt="%f")
        fp.write("\n")
        np.savetxt(fp, data2, fmt="%f")
        fp.close()
      else:
        fp = open(fname, "w")
        fp.write("# 1\n")
        fp.write("# %d"%data1.shape[0])
        for i in range(len(data2)):
          fp.write(":%d"%data2[i].shape[0])
        fp.write("\n")
        np.savetxt(fp, data1, fmt="%f")
        for i in range(len(data2)):
          fp.write("\n")
          np.savetxt(fp, data2[i], fmt="%f")
        
        fp.close()
    else:
      fp = open(fname, "w")
  
      # write header
      fp.write("# %d\n"%len(data1))
      for i in range(len(data1)):
        fp.write("# %d"%data1[i].shape[0])
        if not isinstance(data2[i], list):
          fp.write(":%d\n"%data2[i].shape[0])
        else:
          for j in range(len(data2[i])):
            fp.write(":%d"%data2[i][j].shape[0])
          fp.write("\n")
      # write data
      ic = 0
      for i in range(len(data1)):
        if ic != 0:
          fp.write("\n")
        np.savetxt(fp, data1[i], fmt="%f")
        if not isinstance(data2[i], list):
          fp.write("\n")
          np.savetxt(fp, data2[i], fmt="%f")
        else:
          for j in range(len(data2[i])):
            fp.write("\n")
            np.savetxt(fp, data2[i][j], fmt="%f")
        
        ic += 1
      
      fp.close()
##################################################
def gauss(t, f, tau, wid):
  return f/np.sqrt(2*np.pi*wid**2)*np.exp(-0.5 * (t-tau)**2/(wid)**2)

def tophat(t, f, tau, wid):
  ret = np.zeros(t.shape)
  idx = np.where((t>tau-wid/2) & (t<tau+wid/2))
  ret[idx[0]] = f/wid
  return ret

def gamma(t, f, tau, wid):
  ret = np.zeros(t.shape)
  ret[t>=tau] = f*(t[t>=tau]-tau)/wid**2 * np.exp(-(t[t>=tau]-tau)/wid)
  return ret

def exponential(t, f, tau, wid):
  ret = np.zeros(t.shape)
  ret[t>=tau] = f/wid * np.exp(-(t[t>=tau]-tau)/wid)
  return ret
##################################################

def simlc(tfs={"comp1":("gauss", 1.0, 10.0, 5.0), "comp2":("gauss", 1.0, 30.0, 5.0)},
          lag_range=[-10, 50], doshow=True, tspan=200, dt=2.0):
  """
  generate mock light curves

  """
  print("generate mock light curves...")

  Tot = tspan*2.5
  ts, fs = drw(Tot, [0.3, 50.0, 1.0e-100])
  #ts, fs = psd(Tot, [1.0e-5, 2.5, 1.0e-100])

  fe = np.zeros(ts.shape[0]) + float(0.01)
  con = np.stack((ts, fs, fe), axis=-1)
  dt_org = con[1, 0] - con[0, 0]
  
  ntau = int((lag_range[1]-lag_range[0])/dt_org)
  resp = np.zeros(ntau)
  tau = np.array(np.arange(ntau))*dt_org + lag_range[0]

  funcs = {"gauss":gauss, "tophat": tophat, "gamma": gamma, "exp":exponential}

  for key in tfs.keys():
    tf = tfs[key][0]
    model = tfs[key][1:]
    resp[:] += funcs[tf](tau, model[0], model[1], model[2])

  norm = np.sum(np.abs(resp[:])) * dt_org
  resp[:] /= norm
  conv_org = convolve_fft(con[:, 1]-1.0, resp) * dt_org + 1.0
  # to account for tau[0]!=0
  conv = np.interp(ts, ts+tau[0], conv_org)


  plt.plot(con[:, 0], con[:, 1])
  plt.plot(con[:, 0], conv)
  if doshow:
    plt.show()
  else:
    plt.close()

  np.savetxt("data/resp_input.txt", np.column_stack((tau, resp)))

  DT = dt
  tline0 = 0.0
  tline1 = tspan
  nline = int((tline1-tline0)/DT+0.5)+1
  line = np.zeros((nline, 3))
  line[:, 0] = np.linspace(tline0, tline1, nline)
  line[:, 1] = np.interp(line[:, 0], con[:, 0], conv)
  line[:, 2] = float(0.01)

  # add Gaussian errors
  con[:, 1] += np.random.randn(con.shape[0]) * con[:, 2]
  line[:,1] += np.random.randn(nline) * line[:, 2]

  incr = int(DT/(con[1, 0]-con[0, 0])+0.5)
  if incr < 1:
    incr = 1
  idx = np.where((con[:, 0] >= 0.0) & (con[:, 0] <= tspan))
  con_out = copy.copy(con[idx[0], :])
  con_out = con_out[::incr, :]

  # print(con_out[1, 0]-con_out[0, 0], np.std(con[:, 1]))
  # print(con_out.shape, line.shape)

  fig = plt.figure()

  ax = fig.add_subplot(121)
  plt.errorbar(con_out[:, 0], con_out[:, 1], yerr=con_out[:, 2])
  plt.errorbar(line[:, 0], line[:, 1], yerr = line[:, 2])

  ax = fig.add_subplot(122)
  plt.plot(tau, resp)
  if doshow:
    plt.show()
  else:
    plt.close()

  rnd_con = np.unique(np.random.randint(con_out.shape[0], size=int(con_out.shape[0]*0.8)))
  rnd_line = np.unique(np.random.randint(line.shape[0], size=int(line.shape[0]*0.8)))
  rnd_con = np.sort(rnd_con)
  rnd_line = np.sort(rnd_line)
  format_mica("data/sim.txt", con_out, line)

  print("output light curves to data/sim.txt.")
  print("output transfer function to data/resp_input.txt")


if __name__ == "__main__":
  
  # tfs={"comp1":("gauss", 1.0, -10.0, 5.0), "comp2":("tophat", 1.0, 30.0, 5.0)}
  tfs={"comp1":("gamma", 1.0, 0.0, 6.0), "comp2":("gauss", 0.5, 50.0, 5.0)}
  simlc(tfs, lag_range=[-10, 100], doshow=True, tspan=300, dt=1.0)

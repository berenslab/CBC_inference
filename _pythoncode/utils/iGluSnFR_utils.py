from scipy.signal import convolve
import numpy as np
import warnings

import math_utils
import time_utils


############################################################################
def rate2iGluSnFR(trace, rec_time=None, dt=None, T=0.5, n_drop=0, kind='single'):
  ''' Convolve rate trace with iGluSnFR kernel.
  '''
    
  trace = np.array(trace)
  
  dt = time_utils.get_and_test_dt(dt=dt, t_array=rec_time)
  
  # iGluSnFR kernel.
  if kind == 'double':
    kernel = get_iGluSnFR_kernel_double_exp(duration=T, dt=dt)
  elif kind == 'single':
    kernel = get_iGluSnFR_kernel_exp_decay(duration=T, dt=dt)
  else:
    raise NotImplementedError
  
  # Convolve.
  iGluSnFR = convolve(trace, kernel, mode='full')[0:trace.size]
  
  # Drop values at the beginning were convolution isn't meaningful.
  if n_drop > 0: iGluSnFR[0:n_drop] = iGluSnFR[n_drop]
  
  # Normalize.
  iGluSnFR = math_utils.normalize(iGluSnFR)

  return iGluSnFR
  
############################################################################
def get_iGluSnFR_kernel_double_exp(
      duration, dt,
      tau_r=-0.09919711, # rise time constant
      tau_d=-0.04098927, # decay time constant
    ):  
    ''' Compute double exponential kernel.
    '''
    t_kernel = np.arange(0, duration, dt)
    kernel = np.exp(t_kernel / tau_d) - np.exp((tau_d + tau_r) / (tau_d * tau_r) * t_kernel)
    kernel /= np.max(kernel)
    return kernel.copy()
    
############################################################################
def get_iGluSnFR_kernel_exp_decay(
      duration, dt,
      tau=0.06, # fall time constant
    ):
    ''' Compute exponential decay kernel.
    '''
    kernel = np.exp(-np.arange(0,duration,dt)/0.06)
    return kernel.copy()
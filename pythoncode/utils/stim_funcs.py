import numpy as np
from pandas import DataFrame
from matplotlib import pyplot as plt
from scipy.optimize import minimize

import interpolation_utils

class StimulusGenerator(object):

  ############################################################################
  def __init__(
      self, stim_time, n_params, stim_mode,
      predur=0.0, postdur=0.0, stim_mulitplier=1.0,
      normalize_stim=False, spline_mode='cubic',
      stim_time_min=None,
    ):
    self.stim_mode = stim_mode
    self.stim_time = stim_time
    
    self.stim_time_min = stim_time_min
    
    self.predur = predur
    self.postdur = postdur
    assert stim_time.max() > (predur + postdur)*1.0
    
    # Find indexes in stim time to start and stop stimulus.
    self.idx_start_stim = np.where(stim_time > predur)[0][0]
    post_dur_idxs = np.where((stim_time - stim_time.max()) > -postdur)[0]
    if post_dur_idxs.size == 0: self.idx_stop_stim = self.stim_time.size
    else: self.idx_stop_stim = post_dur_idxs[0]
    
    # Save number of parameters to expect.
    self.n_params = n_params
  
    self.stim_mulitplier = stim_mulitplier
    self.normalize_stim  = normalize_stim
    self.spline_mode     = spline_mode
  
  ############################################################################  
  def create_zero_stimulus(
      self, filename=None, verbose=False, plot=False,
    ):
    ''' Create zero current stimulus.
    Parameters:
      
    filename : str
      Name of output file.
      
    verbose : bool
      Print user information?
      
    plot : bool
      Plot stimulus?
    '''
    stim = np.zeros(self.stim_time.size)
    
    if plot:
      self.plot_stimulus(stim)
    
    if filename is not None:
      self.save_stimulus_to_file(filename, self.stim_time, stim, verbose=verbose)
  
  ############################################################################
  @staticmethod
  def save_stimulus_to_file(filename, stim_time, stim, verbose=True):
    ''' Save given stimulus to given file.
    Will save as csv data.
    '''
    if verbose: print('Save to ' + filename)
    stim_df = DataFrame({'Time': stim_time, 'Stim': stim})
    stim_df.to_csv(filename, index=False, header=False, columns=['Time', 'Stim'])
    
  ############################################################################  
  def create_stimulus(self, params, filename=None, plot=False, verbose=False, **kwargs):
    ''' Create a stimulus.
    '''
    if self.stim_mode == 'spline':
      stim, stim_anchor_points_time, stim_anchor_points_amp =\
        self.create_stimulus_spline(params=params, verbose=verbose, **kwargs)
    
    elif self.stim_mode == 'charge neutral':
      stim, stim_anchor_points_time, stim_anchor_points_amp =\
        self.create_stimulus_charge_neutral(params=params, verbose=verbose, var_dur=False, **kwargs)
    
    elif self.stim_mode == 'charge neutral var dur':
      assert self.stim_time_min is not None, 'Please define.'
      stim, stim_anchor_points_time, stim_anchor_points_amp =\
        self.create_stimulus_charge_neutral(params=params, verbose=verbose, var_dur=True, **kwargs)
    
    else:
      raise NotImplementedError()

    if self.normalize_stim:
      stim_anchor_points_amp /= np.max(np.abs(stim))
      stim /= np.max(np.abs(stim))
    
    stim *= self.stim_mulitplier
    stim_anchor_points_amp *= self.stim_mulitplier

    if plot:
      self.plot_stimulus(stim, stim_anchor_points_time, stim_anchor_points_amp)

    if filename is not None:
      self.save_stimulus_to_file(filename, self.stim_time, stim, verbose=verbose)

    return stim

  ############################################################################  
  def create_stimulus_spline(self, params, verbose=False):
    stim, *_ = self.__create_stim_from_params(params)    
    return stim
  
  ############################################################################  
  def create_stimulus_charge_neutral(
      self, params, filename=None, verbose=False, var_dur=False,
    ):
    
    params = np.asarray(params)
    assert params.size == self.n_params, str(params.size) + '!=' + str(self.n_params)
    
    n_max_iter = 20
    
    if var_dur:
      loss_fun = self.__compute_total_charge_var_dur
    else:
      loss_fun = self.__compute_total_charge

    for i in range(n_max_iter):
      mean_x0 = np.sum(params)
      x_opt = minimize(
        loss_fun,
        x0=np.random.normal(-mean_x0, (0.2*(i+1))*self.n_params), args=(params),
        bounds=[(mean_x0-(i+3)*np.abs(mean_x0)-i, mean_x0+(i+3)*np.abs(mean_x0)+1)], method='SLSQP',
      )
      loss = x_opt.fun
      best_param = x_opt.x
    
      if loss < 1e-4: break
    
    if loss > 1e-4:
      print('WARNING: stimulus is not charge neutral {:.4f}'.format(loss))
    
    stim, stim_anchor_points_time, stim_anchor_points_amp = loss_fun(
      best_param, params, return_charge=False, return_anchors=True,
    )
    
    return stim, stim_anchor_points_time, stim_anchor_points_amp
  
  ############################################################################
  def __compute_total_charge(self, param_x, params, return_charge=True, return_anchors=False):
    
    params = np.asarray(params).flatten()
    assert params.size == self.n_params, str(params.size) + '!=' + str(self.n_params)
    
    stim, stim_anchor_points_time, stim_anchor_points_amp =\
      self.__create_stim_from_params(params=np.append(params, param_x))
    
    if return_charge:
      return np.abs(np.sum(stim))
    elif return_anchors:
      return stim, stim_anchor_points_time, stim_anchor_points_amp
    else:
      return stim
  
  ############################################################################
  def __create_stim_from_params(self, params):

    params = np.asarray(params).flatten()
    assert params.size == self.n_params+1, str(params.size) + '!=' + str(self.n_params)
    
    # Get indexes to create stimulus with parameters. Add two zeros.
    stim_anchor_idxs = np.linspace(self.idx_start_stim, self.idx_stop_stim-1, params.size+2, dtype='int')
    stim_anchor_points_time = self.stim_time[stim_anchor_idxs]
    stim_anchor_points_amp = np.concatenate([np.zeros(1), params, np.zeros(1)])
    
    stim = interpolation_utils.in_ex_polate(
      x_old=stim_anchor_points_time,
      y_old=stim_anchor_points_amp,
      x_new=self.stim_time,
      kind=self.spline_mode
    )
    
    return stim, stim_anchor_points_time, stim_anchor_points_amp
  
  ############################################################################
  def __compute_total_charge_var_dur(self, param_x, params, return_charge=True, return_anchors=False):
  
    params = np.asarray(params).flatten()
    assert params.size == self.n_params
  
    assert (params[-1] >= 0) and (params[-1] <= 1)
  
    stim_time_max = self.stim_time.max() - self.predur - self.postdur
    assert stim_time_max > self.stim_time_min
  
    subtract_stim_time = self.postdur + (1-params[-1]) * (stim_time_max - self.stim_time_min)
  
    post_dur_idxs = np.where((self.stim_time - self.stim_time.max()) > -subtract_stim_time)[0]
    
    if post_dur_idxs.size == 0:
      idx_stop_stim = self.stim_time.size
    else:
      idx_stop_stim = post_dur_idxs[0]
  
    stim_anchor_idxs = np.linspace(self.idx_start_stim, idx_stop_stim-1, params.size+2, dtype='int')
    stim_anchor_points_time = self.stim_time[stim_anchor_idxs]
    stim_anchor_points_amp = np.concatenate([np.zeros(1), params[:-1], param_x*np.ones(1), np.zeros(1)])
    
    stim = interpolation_utils.in_ex_polate(
      x_old=stim_anchor_points_time,
      y_old=stim_anchor_points_amp,
      x_new=self.stim_time,
      kind=self.spline_mode,
    )
    
    if return_charge:
      return np.abs(np.sum(stim))
    elif return_anchors:
      return stim, stim_anchor_points_time, stim_anchor_points_amp
    else:
      return stim

  ############################################################################
  def plot_stimulus(self, stim, stim_anchor_points_time=None, stim_anchor_points_amp=None):
    idx_change = np.argwhere(((stim >= 0).astype(int)[1:] - (stim < 0).astype(int)[0:-1])==0).flatten()
  
    plt.figure(figsize=(12,2))
    plt.plot(self.stim_time, stim)
    for idx in idx_change:
      plt.axvline(self.stim_time[idx], c='k', linestyle='--')
      
    if (stim_anchor_points_time is not None) and (stim_anchor_points_amp is not None):
      for idx, (time, amp) in enumerate(zip(stim_anchor_points_time, stim_anchor_points_amp)):
        if idx == 0 or idx == stim_anchor_points_time.size-1:
          color = 'gray'
        else:
          color = 'red'
        plt.plot(time, amp, marker='x', color=color, ms='10')
      
    plt.axhline(0,c='k')
    plt.title("Charge = {:.3g}".format(np.mean(np.abs(stim)) * (self.stim_time[-1]-self.stim_time[0])))
    plt.ylabel('Stim (uA)')
    plt.xlabel('Time')
    plt.axvline(self.predur, c='r')
    plt.show()
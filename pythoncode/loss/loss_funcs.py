import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import warnings

import interpolation_utils
import lin_trans_utils
import iGluSnFR_utils

import importlib
importlib.reload(lin_trans_utils)

############################################################################
class Loss():
  def __init__(self):
    pass

  def calc_loss(self):
    return None
  
  def print_loss(self, sample_loss):
    print('Loss')
    for loss_name, loss_value in sample_loss.items():
      print('\t' + loss_name.ljust(30) + ' = {:.2g}'.format(loss_value))
    print()

############################################################################
class LossOptimizeCell(Loss):
   
  loss_params_iGluSnFR = {
    'iGluSnFR':               {},
  }
  
  loss_params_CBC_on = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (None, 3),        'acceptable': (None, 7)},
    'rate_rel_range':         {'good': (5, 10),          'acceptable': (0, 20)},
    'Vm_rest':                {'good': (-65e-3, -45e-3), 'acceptable': (-80e-3, -30e-3)},
    'Vm_rel_range':           {'good': (5e-3, 15e-3),    'acceptable': (0e-3, 30e-3)},
    'Vm_low':                 {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'Vm_high':                {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
  }
  
  loss_params_CBC_off = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (None, 4),        'acceptable': (None, 7)},
    'rate_rel_range':         {'good': (5, 10),          'acceptable': (0, 20)},
    'Vm_rest':                {'good': (-65e-3, -45e-3), 'acceptable': (-80e-3, -30e-3)},
    'Vm_rel_range':           {'good': (5e-3, 25e-3),    'acceptable': (0e-3, 40e-3)},
    'Vm_low':                 {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'Vm_high':                {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
  }
  
  loss_params_CBC_on_no_rate_limit = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (None, 3),        'acceptable': (None, 7)},
    'rate_rel_range':         {'good': (5, None),        'acceptable': (0, None)},
    'Vm_rest':                {'good': (-65e-3, -45e-3), 'acceptable': (-80e-3, -30e-3)},
    'Vm_rel_range':           {'good': (5e-3, 15e-3),    'acceptable': (0e-3, 30e-3)},
    'Vm_low':                 {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'Vm_high':                {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
  }
  
  loss_params_CBC_off_rate_limit = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (None, 4),        'acceptable': (None, 7)},
    'rate_rel_range':         {'good': (5, None),        'acceptable': (0, None)},
    'Vm_rest':                {'good': (-65e-3, -45e-3), 'acceptable': (-80e-3, -30e-3)},
    'Vm_rel_range':           {'good': (5e-3, 25e-3),    'acceptable': (0e-3, 40e-3)},
    'Vm_low':                 {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'Vm_high':                {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
  }
  
  loss_params_cone = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (50, 80),         'acceptable': (0, 100)},
    'rate_rel_range':         {'good': (50, 65),         'acceptable': (0, 100)},
    'Vm_rest':                {'good': (-55e-3, -40e-3), 'acceptable': (-80e-3, -20e-3)},
    'Vm_rel_range':           {'good': (5e-3, 10e-3),    'acceptable': (0e-3, 20e-3)},
    'Vm_low':                 {'good': (-60e-3, None),   'acceptable': (-75e-3, None)},
    'Vm_high':                {'good': (None, -35e-3),   'acceptable': (None, -20e-3)},
  }
  
  
  loss_params_ac = {
    'iGluSnFR':               {},
    'rate_rest':              {'good': (None, 4),        'acceptable': (None, 7)},
    'rate_rel_range':         {'good': (5, None),        'acceptable': (0, None)},
    'Vm_rest':                {'good': (-65e-3, -45e-3), 'acceptable': (-80e-3, -30e-3)},
    'Vm_rel_range':           {'good': (5e-3, 25e-3),    'acceptable': (0e-3, 40e-3)},
    'Vm_low':                 {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'Vm_high':                {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
    'AC_Vm_low':              {'good': (-80e-3, None),   'acceptable': (-100e-3, None)},
    'AC_Vm_high':             {'good': (None, -10e-3),   'acceptable': (None, 0e-3)},
  }

  
  ############################################################################
  def __init__(
      self, rec_time, target, loss_params, t_drop=0,
      mode='lin', absolute=True
    ):

    self.select_loss_params(loss_params) 

    self.rec_time = np.asarray(rec_time)   

    self.t_drop = t_drop
    assert self.t_drop < self.rec_time[-1], 't_drop greater than rec_time.'

    self.init_target(target)
    self.mode = mode
    self.absolute = absolute
  
  ############################################################################   
  def select_loss_params(self, loss_params):
    if loss_params == 'BC ON':
      self.loss_params = self.loss_params_CBC_on
    elif loss_params == 'iGlu only':
      self.loss_params = self.loss_params_iGluSnFR
    elif loss_params == 'BC OFF':
      self.loss_params = self.loss_params_CBC_off
    elif loss_params == 'BC ON no rate limit':
      self.loss_params = self.loss_params_CBC_on_no_rate_limit
    elif loss_params == 'BC OFF no rate limit':
      self.loss_params = self.loss_params_CBC_off_rate_limit
    elif loss_params == 'AC':
      self.loss_params = self.loss_params_ac
    elif loss_params == 'Cone eq and range':
      self.loss_params = self.loss_params_cone
    else:
      raise NotImplementedError()
  
  ############################################################################   
  def init_target(self, target):
    ''' Initialize target, i.e. save target and interpolate for rec_time.
    '''
    assert isinstance(target, pd.DataFrame)
    assert 'Time' and 'mean' in target.columns
    
    self.target_original = target.copy()
    self.update_target()
  
  ############################################################################
  def update_target(self, rec_time=None):
    ''' Interpolate target for given recording time.
    '''
    
    if rec_time is not None:
      self.rec_time = np.asarray(rec_time)
    
    self.n_drop_trace = self.get_n_drop(self.rec_time, self.t_drop)
    
    orignal_time = self.target_original['Time'].values
    self.target_time = orignal_time[(orignal_time >= self.rec_time[0]) & (orignal_time <= self.rec_time[-1])]
    
    self.n_drop_target = self.get_n_drop(self.target_time, self.t_drop)
    
    target_release = self.target_original['mean'].values
    
    self.target = interpolation_utils.in_ex_polate(
      x_old=orignal_time, y_old=target_release, x_new=self.target_time,
    )
  
  ############################################################################
  @staticmethod
  def get_n_drop(time, t_drop):
    return np.argmax(time >= t_drop)
  
  ############################################################################
  def calc_loss(self, rec_data, verbose=False, plot=False):
    
    if rec_data is None:
      sample_loss = {loss_name: np.nan for loss_name in self.loss_params}
      sample_loss['total'] = np.nan
      if verbose: self.print_loss(sample_loss)
      return sample_loss
    
    # Compute loss.
    sample_loss = {}
    for name in self.loss_params:

      if   name == 'iGluSnFR':
        sample_loss[name] = self.loss_iGluSnFR(rec_data['rate'])        
      else:
        
        if name == 'rate_rest':
          value=rec_data['rate'][self.n_drop_trace]
        elif name == 'rate_rel_range':
          value=np.max(np.abs(rec_data['rate'][self.n_drop_trace:] - rec_data['rate'][self.n_drop_trace]))
        elif name == 'Vm_rest':
          value=rec_data['Vm'][self.n_drop_trace]
        elif name == 'Vm_rel_range':
          value=np.max(np.abs(rec_data['Vm'][self.n_drop_trace:] - rec_data['Vm'][self.n_drop_trace]))
        elif name == 'Vm_low':
          value=np.min(rec_data['Vm'][self.n_drop_trace:])                                     
        elif name == 'Vm_high':
          value=np.max(rec_data['Vm'][self.n_drop_trace:])     
        elif name == 'AC_Vm_low':
          value=np.min(rec_data['Vm AC'][self.n_drop_trace:])
        elif name == 'AC_Vm_high':
          value=np.max(rec_data['Vm AC'][self.n_drop_trace:])                
        else:
          raise NotImplementedError()
          
        sample_loss[name] = self.loss_value_in_range(
          value=value,
          good=self.loss_params[name]['good'],
          acceptable=self.loss_params[name]['acceptable']
        )
    
    sample_loss['total'] = np.sum(np.abs(list(sample_loss.values())))
    
    if verbose: self.print_loss(sample_loss)
    if plot: self.plot_sample_loss(rec_data, sample_loss)
    
    return sample_loss
  
  ############################################################################
  def loss_iGluSnFR(self, trace):
    assert self.rec_time.size == trace.size, str(self.rec_time.size) + '!=' + str(trace.size)
  
    shifted_iGluSnFR_trace, iGluSnFR_loss = self.rate2best_iGluSnFR_trace(trace=trace)
    _, f_norm_loss = self.rate2best_iGluSnFR_trace(trace=np.zeros(trace.size))
    
    iGluSnFR_loss /= f_norm_loss
  
    if iGluSnFR_loss > 1:
      if not np.isclose(iGluSnFR_loss, 1):
        plt.figure(figsize=(12,3))
        plt.plot(shifted_iGluSnFR_trace)
        plt.show()
        print("Should be <= 1 by definition but is", iGluSnFR_loss)
      iGluSnFR_loss = 1.0
      
    return iGluSnFR_loss
  
  ############################################################################
  def rate2best_iGluSnFR_trace(self, trace):
  
    iGluSnFR_trace = iGluSnFR_utils.rate2iGluSnFR(
      trace, rec_time=self.rec_time, n_drop=self.n_drop_trace
    )
    
    intpol_iGluSnFR_trace = interpolation_utils.in_ex_polate(
      x_old=self.rec_time, y_old=iGluSnFR_trace, x_new=self.target_time
    )
    
    trans_iGluSnFR_trace, iGluSnFR_loss = lin_trans_utils.best_lin_trans(
      trace=intpol_iGluSnFR_trace, target=self.target, loss_fun=self.compute_iGluSnFR_trace_loss
    )
    
    return trans_iGluSnFR_trace, iGluSnFR_loss
  
  ############################################################################
  def compute_iGluSnFR_trace_loss(self, trace, target):
    return np.mean((trace[self.n_drop_target:] - target[self.n_drop_target:])**2)
  
  ############################################################################
  # Penalize release rates greater then upper bound.
  def loss_value_in_range(self, value, good, acceptable):
    if self.mode == 'linear':
      loss = self.loss_value_in_range_lin(value, good, acceptable)
    elif self.mode == 'sinus':
      loss = self.loss_value_in_range_sin(value, good, acceptable)
    elif self.mode == 'gauss':
      loss = self.loss_value_in_range_gauss(value, good, acceptable)
    else:
      raise NotImplementedError(self.mode)
    
    if self.absolute: loss = np.abs(loss)
    return loss
  
  ############################################################################
  # Penalize release rates greater then upper bound.
  @staticmethod
  def loss_value_in_range_lin(value, good, acceptable):
  
    # Check if too small.
    if good[0] is not None:
      if value <= acceptable[0]:
        loss = -1.0
        return loss
      elif value < good[0]:
        loss = (value - good[0]) / (good[0] - acceptable[0])
        return loss
      else:
        pass
    
    # Check if too large.
    if good[1] is not None:
      if value >= acceptable[1]:
        loss = 1.0
        return loss
      elif value > good[1]:
        loss = (value - good[1]) / (acceptable[1] - good[1])
        return loss
      else:
        pass

    return 0.0
    
  ############################################################################
  # Penalize release rates greater then upper bound.
  @staticmethod
  def loss_value_in_range_sin(value, good, acceptable):
    
    # Check if too small.
    if good[0] is not None:
      if value < good[0]:
        half_point = good[0] - (good[0]-acceptable[0])/2
        if value > half_point:
          loss = -0.5+0.5*np.cos((value - good[0]) / (good[0] - acceptable[0]) * np.pi)
          return loss
        else:
          loss = -1+1/(1+np.exp((value-half_point) / (acceptable[0] - good[0])* 2*np.pi))
          return loss
      else:
        pass # Looks good.
    
    # Check if too large.
    if good[1] is not None:
      
      if value > good[1]:
        half_point = good[1] + (acceptable[1]-good[1])/2
        if value < half_point:
          loss = 0.5-0.5*np.cos((value - good[1]) / (acceptable[1] - good[1]) * np.pi)
          return loss
        else:
          loss = 1-1/(1+np.exp((value - half_point) / (acceptable[1] - good[1])* 2*np.pi))
          return loss
      else:
        pass # Looks good. 
    
    return 0.0
    
  ############################################################################
  # Penalize release rates greater then upper bound.
  @staticmethod
  def loss_value_in_range_gauss(value, good, acceptable):
    
    # Check if too small.
    if good[0] is not None:
      if value < good[0]:
        sigma = (good[0]-acceptable[0])/2
        loss = -1 + np.exp( -(value-good[0])**2 / (2*sigma**2) ) 
        return loss
      else:
        pass # Looks good.
    
    # Check if too large.
    if good[1] is not None:
      
      if value > good[1]:
        sigma = (acceptable[1]-good[1])/2
        loss = 1 - np.exp( -(good[1]-value)**2 / (2*sigma**2) ) 
        return loss
      else:
        pass # Looks good. 
    
    return 0.0
        
  ############################################################################
  def plot_loss_params(self):
        
    fig, axs = plt.subplots(len(self.loss_params),1,figsize=(12,len(self.loss_params)*1.5))

    for ax, (loss_name, loss_dict) in zip(axs, self.loss_params.items()):
      ax.set_title(loss_name)
    
      if ('good' in loss_dict) and ('acceptable' in loss_dict):
        
        xticks = []
        lb = None
        ub = None
      
        if loss_dict['good'][0] is not None:
          lb_g = loss_dict['good'][0]
          xticks.append(lb_g)
        if loss_dict['acceptable'][0] is not None:
          lb_a = loss_dict['acceptable'][0]
          xticks.append(lb_a)            
        if (loss_dict['good'][0] is not None) and (loss_dict['acceptable'][0] is not None):
          lb_rng = np.abs(lb_g-lb_a)
          lb = lb_a - 0.5*lb_rng
          
        if loss_dict['good'][1] is not None:
          ub_g = loss_dict['good'][1]
          xticks.append(ub_g)
        if loss_dict['acceptable'][1] is not None:
          ub_a = loss_dict['acceptable'][1]
          xticks.append(ub_a)
        if (loss_dict['good'][1] is not None) and (loss_dict['acceptable'][1] is not None):
          ub_rng = np.abs(ub_g-ub_a)
          ub = ub_a + 0.5*ub_rng
      
        if lb is None:
          lb = ub - 2*ub_rng
        elif ub is None:
          ub = lb + 2*lb_rng
      
        in_values = np.linspace(lb, ub, 100)
        out_values = np.full(in_values.size, np.nan)
      
        for idx, value in enumerate(in_values):
          out_values[idx] = self.loss_value_in_range(
          value=value, good=loss_dict['good'], acceptable=loss_dict['acceptable']
        )
        
        ax.plot(in_values, out_values)
        ax.set_xticks(xticks)
        
        for xtick in xticks:
          ax.axvline(xtick, c='gray', ls='--')
        if not self.absolute:
          ax.axhline(0, c='gray', ls='--')
        
      elif 'iGluSnFR' in loss_name:
        plot_losses_dict = {}
        plot_losses_dict['sinus'] = self.loss_iGluSnFR(np.sin(10*self.rec_time))
        plot_losses_dict['zeros'] = self.loss_iGluSnFR(np.zeros(self.rec_time.size))
        plot_losses_dict['ones']  = self.loss_iGluSnFR(np.zeros(self.rec_time.size))
        plot_losses_dict['slope'] = self.loss_iGluSnFR(-np.arange(self.rec_time.size))
        plot_losses_dict['noise'] = self.loss_iGluSnFR(np.random.uniform(-1,1,self.rec_time.size))
        
        interpol_target = interpolation_utils.in_ex_polate(self.target_time, self.target, self.rec_time)
        
        plot_losses_dict['noisy target'] = self.loss_iGluSnFR(interpol_target+np.random.normal(0,np.std(self.target),self.rec_time.size))
        plot_losses_dict['flipped target'] = self.loss_iGluSnFR(-interpol_target)

        ax.set_yscale('log')
        ax.bar(np.arange(len(plot_losses_dict)), plot_losses_dict.values())
        ax.set_xticks(np.arange(len(plot_losses_dict)))
        ax.set_xticklabels(plot_losses_dict.keys())
        for idx, plot_losses_value in enumerate(plot_losses_dict.values()):
          ax.text(idx, 1, "{:.4f}".format(plot_losses_value), ha='center', va='top')
            
    plt.tight_layout()
    plt.show()
        
  ############################################################################
  def plot_sample_loss(self, rec_data, sample_loss):
    
    plot_iGluSnFR = 'iGluSnFR' in self.loss_params
    
    n_subplots = 2 + 2*plot_iGluSnFR
    
    fig, axs = plt.subplots(n_subplots, 1, figsize=(12,6), sharex=True)
    
    # Plot rate.
    axs[0].set_title('rate')
    axs[0].plot(self.rec_time, rec_data['rate'])
    axs[0].axvline(self.rec_time[self.n_drop_trace], color='k', alpha=0.3)
    
    yticks = [rec_data['rate'].min(), rec_data['rate'][self.n_drop_trace], rec_data['rate'].max()]
    axs[0].set_yticks(yticks)
    axs[0].set_yticklabels(["{:.2g}".format(ytick) for ytick in yticks])
    for ytick in yticks:
      axs[0].axhline(ytick, color='k', alpha=0.3)
    axs[0].set_ylabel('Rate (ves/s)')
    
    # Plot Vm.
    axs[1].set_title('Vm')
    axs[1].plot(self.rec_time, 1000*rec_data['Vm'])
    axs[1].axvline(self.rec_time[self.n_drop_trace], color='k', alpha=0.3)
    
    yticks = [1000*rec_data['Vm'].min(), 1000*rec_data['Vm'][self.n_drop_trace], 1000*rec_data['Vm'].max()]
    axs[1].set_yticks(yticks)
    axs[1].set_yticklabels(["{:.2g}".format(ytick) for ytick in yticks])
    for ytick in yticks:
      axs[1].axhline(ytick, color='k', alpha=0.3)
    axs[1].set_ylabel('Voltage (mV)')
    
    # Plot iGluSnFR.
    if plot_iGluSnFR:
      axs[2].set_title('iGluSnFR')
      axs[2].plot(self.target_time, self.rate2best_iGluSnFR_trace(rec_data['rate'])[0], c='r', label='shifted trace')
      axs[2].plot(self.target_time, self.target, c='k', alpha=0.5, label='target')
      axs[2].axvline(self.rec_time[self.n_drop_trace], color='k', alpha=0.3)
      axs[2].legend(loc='upper right')
      
      axs[3].set_title(r'$\Delta$ iGluSnFR')
      difference = self.rate2best_iGluSnFR_trace(rec_data['rate'])[0] - self.target
      axs[3].fill_between(self.target_time[self.n_drop_target:], difference[self.n_drop_target:], color='r', label='included', alpha=0.6)
      axs[3].fill_between(self.target_time[:self.n_drop_target], difference[:self.n_drop_target], color='k', label='exluded', alpha=0.3)
      axs[3].axvline(self.target_time[self.n_drop_target], color='k', alpha=0.3)
      axs[3].legend(loc='upper right')
    
    axs[-1].set_xlabel('Time')
    
    plt.tight_layout()
    plt.show()
        
############################################################################
############################################################################

class LossOptimizeAC(LossOptimizeCell):
  loss_params_optimize_ac = {
    'iGluSnFR':           {'max': 1, 'w': 1},
    'max_release_rate':   {'max': 1, 'w': 1, 'upper': 600},
    'release_rate_range': {'max': 1, 'w': 1, 'lower_absolut': 4,     'lower_relative': 0.25},
    'release_rate_mean':  {'max': 1, 'w': 1, 'lower': 2,             'upper': 12},
    
    'Vm_range_AC':        {'max': 1, 'w': 1, 'lower_absolut': 10e-3, 'lower_relative': 0.0},
    'Vm_range_OFF':       {'max': 1, 'w': 1, 'lower_absolut': 10e-3, 'lower_relative': 0.0},
    
    'Vm_limits_AC':       {'max': 1, 'w': 1, 'lower': -100e-3,       'upper': -10e-3},
    'Vm_limits_OFF':      {'max': 1, 'w': 1, 'lower': -100e-3,       'upper': -10e-3},
    
    'Vm_rest_OFF':        {'max': 1, 'w': 1, 'lower': -65e-3,        'upper': -52e-3},
  }
  
  ############################################################################
  def calc_loss(self, rec_data, verbose=False, plot=False):
    
    if rec_data is None:
      sample_loss = {loss_name: np.nan for loss_name in self.loss_params}
      sample_loss['total'] = np.nan
      return sample_loss
    
    # Compute loss.
    sample_loss = {}
    for loss_name in self.loss_params:
      if   loss_name == 'iGluSnFR':           sample_loss[loss_name] = self.loss_iGluSnFR(rec_data['rate OFF'])
      elif loss_name == 'max_release_rate':   sample_loss[loss_name] = self.loss_max_release_rate(rec_data['rate OFF'][self.n_drop_trace:].max())
      elif loss_name == 'release_rate_range': sample_loss[loss_name] = self.loss_release_rate_range(rec_data['rate OFF'][self.n_drop_trace:].min(), rec_data['rate OFF'][self.n_drop_trace:].max())
      elif loss_name == 'release_rate_mean':  sample_loss[loss_name] = self.loss_release_rate(rec_data['rate OFF'][self.n_drop_trace:].mean())
      elif loss_name == 'Vm_range_AC':        sample_loss[loss_name] = self.loss_Vm_range(rec_data['Vm AC'][self.n_drop_trace:].min(), rec_data['Vm AC'][self.n_drop_trace:].max(), name='Vm_range_AC')
      elif loss_name == 'Vm_range_OFF':       sample_loss[loss_name] = self.loss_Vm_range(rec_data['Vm OFF'][self.n_drop_trace:].min(), rec_data['Vm OFF'][self.n_drop_trace:].max(), name='Vm_range_OFF')
      elif loss_name == 'Vm_limits_AC':       sample_loss[loss_name] = self.loss_Vm_limits(rec_data['Vm AC'][self.n_drop_trace:].min(), rec_data['Vm AC'][self.n_drop_trace:].max(), name='Vm_limits_AC')
      elif loss_name == 'Vm_limits_OFF':      sample_loss[loss_name] = self.loss_Vm_limits(rec_data['Vm OFF'][self.n_drop_trace:].min(), rec_data['Vm OFF'][self.n_drop_trace:].max(), name='Vm_limits_OFF')
      elif loss_name == 'Vm_rest_OFF':        sample_loss[loss_name] = self.loss_Vm_rest(rec_data['Vm OFF'][self.n_drop_trace], name='Vm_rest_OFF')
    
    sample_loss['total'] = sum(sample_loss.values())
    
    # Print?
    if verbose: self.print_loss(sample_loss)

    # Plot?
    if plot: self.__plot_loss(rec_data, sample_loss)
    
    return sample_loss
    
  ############################################################################
  def __plot_loss(self, rec_data, sample_loss):
    plt.figure(figsize=(12,6))
    plt.subplot(311)
    plt.title('rate OFF')
    plt.plot(self.rec_time, rec_data['rate OFF'])
    if 'release_rate_mean' in self.loss_params:
      plt.axhline(self.loss_params['release_rate_mean']['lower'], c='r', alpha=0.3)
      plt.axhline(self.loss_params['release_rate_mean']['upper'], c='r', alpha=0.3)
    plt.axvline(self.rec_time[self.n_drop_trace])
    
    plt.subplot(312)
    plt.title('Vm AC')
    plt.plot(self.rec_time, rec_data['Vm AC'])
    plt.axvline(self.rec_time[self.n_drop_trace])
    
    if 'iGluSnFR' in self.loss_params:
      plt.subplot(313)
      plt.title('iGluSnFR')
      plt.plot(self.target_time, self.rate2best_iGluSnFR_trace(rec_data['rate OFF'])[0], 'r', label='target')
      plt.plot(self.target_time, self.target, 'k')
      plt.axvline(self.rec_time[self.n_drop_trace])
    
    plt.xlabel('Time')
    
    plt.tight_layout()
    plt.show()
    
############################################################################
############################################################################
class LossOptimizeCOMSOL_EDL(Loss):
    
  ############################################################################
  def calc_loss(self, rec_data, plot=False, verbose=False):
  
    sample_loss = {}
  
    for f, rec_data_f in rec_data.items():
      for I, rec_data_f_I in rec_data_f.items():
        rec_time = rec_data_f_I['Time'].values
        trace    = rec_data_f_I['Voltage'].values

        sample_loss['f='+str(f)+'_'+str(I)] = I/1000. - np.max(trace)
        
        # Plot?
        if plot: self.__plot_loss(self, rec_time, trace, f, I)
    
    sample_loss['total'] = np.sum(list(sample_loss.values()))
    
    # Print?
    if verbose: self.print_loss(sample_loss)
    
    return sample_loss
  
  ############################################################################
  @staticmethod
  def __plot_loss(rec_time, trace, f, I):
    plt.figure(figsize=(12,1))
    plt.title('f='+str(f)+'_'+str(I))
    plt.plot(rec_time, trace, label='trace')
    plt.legend()
    plt.show()
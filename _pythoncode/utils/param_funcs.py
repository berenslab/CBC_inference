import numpy as np
import warnings
from matplotlib import pyplot as plt

############################################################################
class DummyParameters():
    def __init__(self, n_params_stim):
        self.p_N = n_params_stim
        self.p_names = ['p_'+str(i) for i in range(self.p_N)]
        self.p_use_log = {}
        
    def opt_param2sim_param(self, x, name):
        return x
    
    def sim_param2opt_param(self, x, name):
        return x

############################################################################
class Parameters():

  ############################################################################
  def __init__(
    self,
    p_names          = None,
    p_default        = {},
    p_use_log        = [],
    p_lower_is_nan   = {},
    p_larger_is_nan  = {},
    p_lower_clip     = {},
    p_larger_clip    = {},
    p_range          = {},
    positive_prefix  = None,
  ):
  
    # Use default parameters are param names.
    if p_names is None:
      if len(p_default) > 0:
        p_names = list(p_default.keys())
  
    assert isinstance(p_names, list)
    assert isinstance(p_default, dict)
    assert isinstance(p_range, dict)
    assert isinstance(p_use_log, list)
    assert isinstance(p_lower_is_nan, dict)
    assert isinstance(p_larger_is_nan, dict)
    assert isinstance(p_lower_clip, dict)
    assert isinstance(p_larger_clip, dict)

    # Set optimization parameters.
    self.p_names         = p_names
    self.p_default       = p_default
    self.p_range         = p_range
    self.p_larger_is_nan = p_larger_is_nan
    self.p_lower_is_nan  = p_lower_is_nan
    self.p_lower_clip    = p_lower_clip
    self.p_larger_clip   = p_larger_clip
    self.p_use_log       = p_use_log
    
    self.p_N = len(self.p_names)
    
    # Check if default params are in p_names in within range.
    
    for p_name in p_default.keys():
      assert p_name in p_names, p_name
      if p_name in p_range.keys():
        assert p_default[p_name] >= p_range[p_name][0], p_name
        assert p_default[p_name] <= p_range[p_name][1], p_name
      else:
        assert p_name in p_use_log, p_name
    
    for p_name in p_range.keys():
      assert p_name in p_default.keys(), p_name 
      
    if positive_prefix is not None:
      self.enforce_positive_prefix(positive_prefix)
     
  ############################################################################
  def enforce_positive_prefix(self, suffix):
    for param in self.p_names:
      if param[:len(suffix)] == suffix and (param not in self.p_lower_clip.keys()):
        self.p_lower_clip[param] = 0.0
  
  ############################################################################
  def opt_param2sim_param(self, opt_param, name=None, verbose=False):

    if isinstance(opt_param, np.ndarray):
      opt_param = opt_param.copy().astype(float)
    else:
      opt_param = np.atleast_1d(float(opt_param))
  
    # Find out what to do with parameter.
    if name is None:
      warnings.warn('No parameter name given. Will use default (=no) transformation!')
      use_log = False
      use_range = False
    else:
      assert name in self.p_names
      use_log   = name in self.p_use_log
      use_range = name in self.p_range
    assert not(use_log) or not(use_range), 'Can not use log space and range.'
    
    # Get new value.
    if use_log:
      try:
        sim_param = 2**opt_param
      except:
        if verbose: print('Could not compute power. Return NaN.')
        return np.nan
    elif use_range:
      sim_param = (opt_param * (self.p_range[name][1] - self.p_range[name][0])) + self.p_range[name][0]
    else:
      sim_param = opt_param

    # Check other conditions.
    if name in self.p_larger_is_nan:
      sim_param[sim_param >= self.p_larger_is_nan[name]] = np.nan
    
    if name in self.p_lower_is_nan:
      sim_param[sim_param <= self.p_lower_is_nan[name]] = np.nan
    
    if name in self.p_lower_clip:
      sim_param[sim_param <= self.p_lower_clip[name]] = self.p_lower_clip[name]
      
    if name in self.p_larger_clip:
      sim_param[sim_param >= self.p_larger_clip[name]] = self.p_larger_clip[name]
    
    if sim_param.size == 1:
      sim_param = float(sim_param)
    
    return sim_param
  
  ############################################################################
  def opt_params2sim_params(self, opt_params):
    if opt_params is None: return {}
    opt_params = np.asarray(opt_params)
    assert opt_params.size == len(self.p_names), 'Error: ' + str(opt_params.size) + ' != ' + str(len(self.p_names)) + ' ' + str(self.p_names)
    return {p_name: self.opt_param2sim_param(opt_params[p_idx], name=p_name) for p_idx, p_name in enumerate(self.p_names)}
  
  ############################################################################
  def sim_param2opt_param(self, sim_param, name=None, verbose=False):
    if isinstance(sim_param, np.ndarray):
      sim_param = sim_param.copy().astype(float)
    else:                                 
      sim_param = np.atleast_1d(float(sim_param))
  
    # Find out what to do with parameter.
    if name is None:
      warnings.warn('No parameter name given. Will use default (=no) transformation!')
      use_log = False
      use_range = False
    else:
      use_log = name in self.p_use_log
      use_range = name in self.p_range
      assert not(use_log) or not(use_range), 'Can not use log space and range.'
  
    # Check other conditions.
    if name in self.p_larger_is_nan:
      sim_param[sim_param >= self.p_larger_is_nan[name]] = np.nan
    
    if name in self.p_lower_is_nan:
      sim_param[sim_param <= self.p_lower_is_nan[name]] = np.nan
    
    if name in self.p_lower_clip:
      sim_param[sim_param <= self.p_lower_clip[name]] = 0.0
      
    if name in self.p_larger_clip:
      sim_param[sim_param >= self.p_larger_clip[name]] = 1.0
    
    # Get new value.
    if use_log:
      if verbose: print('use log.')
      is_zero = (sim_param == 0)
      opt_param = sim_param
      if np.any(is_zero):
        opt_param[is_zero] = np.log2(self.p_lower_clip[name])
      opt_param[~is_zero] = np.log2(sim_param)[~is_zero]
    elif use_range:
      if verbose: print('use range.')
      opt_param = (sim_param - self.p_range[name][0]) / (self.p_range[name][1] - self.p_range[name][0])
    else:
      if verbose: print('use none.')
      opt_param = sim_param
      
    if opt_param.size == 1:
      opt_param = float(opt_param)
    
    return opt_param
    
  ############################################################################
  def sim_params2opt_params(self, sim_params):
    if sim_params is None: return []
    assert len(sim_params) == len(self.p_names), 'Error: ' + str(len(sim_params)) + ' != ' + str(len(self.p_names)) + ' ' + str(self.p_names)
    return [self.sim_param2opt_param(sim_params[p_name], name=p_name, verbose=False) for p_idx, p_name in enumerate(self.p_names)]
    
    
  ############################################################################
  def plot(self, opt_bounds=None):
    
    if self.p_N <= 6:
      sp_nx = self.p_N
      sp_ny = 1
    else:
      sp_nx = int(np.ceil(np.sqrt(self.p_N)))
      sp_ny = sp_nx
    
    plt.figure(figsize=(sp_nx*2, sp_ny*2))
    for p_idx, p_name in enumerate(self.p_names):
      
      if opt_bounds is None:
        # Get lower bound.
        if p_name in self.p_range:
          lower = self.p_range[p_name][0]
        elif p_name in self.p_lower_is_nan:
          lower = self.p_lower_is_nan[p_name]
        elif p_name in self.p_lower_clip:
          lower = self.p_lower_clip[p_name]
        elif p_name in self.p_use_log:
          lower = 1e-6
        else:
          lower = 0.0
        
        # Get upper bound.
        if p_name in self.p_range:
          upper = self.p_range[p_name][1]
        elif p_name in self.p_larger_is_nan:
          upper = self.p_larger_is_nan[p_name]
        else:
          upper = 1.0
          
        sim_param_values = np.linspace(lower, upper, 100, endpoint=True)
        opt_param_values = self.sim_param2opt_param(sim_param_values, name=p_name)
      else:
        opt_param_values = np.linspace(opt_bounds[0], opt_bounds[1], 100, endpoint=True)
        sim_param_values = self.opt_param2sim_param(opt_param_values, name=p_name)

      # Get default.
      has_default = p_name in self.p_default
      
      if has_default:
        sim_param_default = self.p_default[p_name]
        opt_param_default = self.sim_param2opt_param(sim_param_default, name=p_name)
      
      plt.subplot(sp_ny, sp_nx, p_idx+1)
      plt.plot(opt_param_values, sim_param_values)
      plt.xticks([opt_param_values[0], opt_param_values[-1]])
      plt.yticks([sim_param_values[0], sim_param_values[-1]])
      
      if has_default:
        plt.plot(opt_param_default, sim_param_default, '*r')
      plt.xlabel('opt space')
      plt.ylabel('sim space')
      plt.title(p_name)
      
    plt.tight_layout()
    plt.show()
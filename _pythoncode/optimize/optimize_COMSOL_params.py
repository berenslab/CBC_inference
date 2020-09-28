import os
from time import sleep
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import warnings

from optim_funcs import Optimizer
import data_utils
import interpolation_utils

############################################################################
############################################################################
class OptimizerCOMSOLparams(Optimizer): 

  def __init__(self, params, output_folder, EDL_phase_total, V_amps, absZ_est, verbose=True, p_unit={}, reset=True):
    assert isinstance(p_unit, dict)
  
    self.output_folder = output_folder
    self.params = params
    self.p_unit = p_unit
    self.EDL_phase_total = EDL_phase_total
    self.V_amps = V_amps
    self.absZ_est = absZ_est
  
    data_utils.make_dir('optim_data/' + output_folder)
    data_utils.make_dir('COMSOL_input')
    data_utils.make_dir('COMSOL_output')
 
    if reset: self.reset_and_wait_for_COMSOL()
  
  ############################################################################
  @staticmethod
  def clear(input=True, output=True):
    if input:
      COMSOL_inputs = [COMSOL_input for COMSOL_input in os.listdir("COMSOL_input")]
      for COMSOL_input in COMSOL_inputs: os.remove("COMSOL_input/" + COMSOL_input)
   
    if output:
      COMSOL_outputs = [COMSOL_output for COMSOL_output in os.listdir("COMSOL_output")]
      for COMSOL_output in COMSOL_outputs: os.remove("COMSOL_output/" + COMSOL_output)
  
  ############################################################################
  @staticmethod
  def set_input_ready():
    with open('_ready_input', 'w') as f:
      f.write(str(1))
  
  ############################################################################
  @staticmethod
  def set_output_ready():
    with open('_ready_output', 'w') as f:
      f.write(str(1))
      
  ############################################################################
  @staticmethod
  def set_input_not_ready():
    if os.path.isfile('_ready_input'):
      os.remove('_ready_input')
      
  ############################################################################
  @staticmethod
  def set_output_not_ready():
    if os.path.isfile('_ready_output'):
      os.remove('_ready_output')
      
  ############################################################################
  @staticmethod
  def is_output_ready():
    if '_ready_output' in os.listdir('.'):
      return True
    else:
      return False
  
  ############################################################################
  def reset_and_wait_for_COMSOL(self):
    self.clear()
    self.set_output_ready()
    print('Please (re-)start local notebook for COMSOL ... ', end='')
    while (self.is_output_ready()):
      sleep(0.3)
      self.set_input_not_ready()
      
    sleep(1)
    self.set_input_not_ready()
      
    print('Ready!')
  
  ############################################################################
  def run_COMSOL(self, verbose=True):
    self.set_input_ready()
    if verbose: print('Wait Vext', end=', ')
    while (not self.is_output_ready()):
      sleep(0.3)
    self.set_output_not_ready()
    if verbose: print('Ready!')
  
  ############################################################################
  def get_unit(self, param):
    if param in self.p_unit:
      return self.p_unit[param]
    else:
      return 1.0
  
  ############################################################################
  def RC(self, phi, f):
    return -np.tan(phi/180*np.pi)/self._w(f)
  
  ############################################################################
  def R_from_CR(self, absZ, RC, f):
    return absZ * np.sqrt(1 + self._w(f)**2 * RC**2)
  
  ############################################################################
  def _w(self, f):
    assert f in [25, 40]
    return (2*np.pi*f)
  
  ############################################################################
  def read_outputs(self, verbose=False):
  
    rec_data = {}
    for f0 in fs:
      rec_data[f0] = {}
      for V0 in enumerate(self.V_amps[f0]):
        filename = 'COMSOL_output/v_EDL_RC_w_f0_' + str(f0) + '_V0_' + str(V0/1000) + '.txt' 
        if verbose: print('Read: ' + filename)
        try:
          rec_data[f0][V0] = pd.read_csv(filename,  comment='%', header=None, names=['Time', 'Current'], delim_whitespace=True)
        except:
          warnings.warn('Could not read ' + filename)
          rec_data[f0][V0] = None
        
    return rec_data
  
  ############################################################################
  def create_RC_inputs(self):
  
    n_rows = len(self.V_amps[25]) + len(self.V_amps[40])
    n_cols = 3
    
    # Prepare. First column = f0 [Hz], second column = V0 [mV]
    Rs = np.zeros((n_rows, n_cols))
    
    Rs[0:len(self.V_amps[25]),0] = 25
    Rs[len(self.V_amps[25]):,0]  = 40
    Rs[0:len(self.V_amps[25]),1] = self.V_amps[25]
    Rs[len(self.V_amps[25]):,1]  = self.V_amps[40]
    
    Cs = Rs.copy()
    
    # Fill arrays.
    for if0, f0 in enumerate([25, 40]):
      for iV0, V0 in enumerate(self.V_amps[f0]):
        # Compute RC through phase.
        RC0 = self.RC(self.EDL_phase_total[f0][V0], f0)
    
        # Get R from RC and |Z|
        R = self.R_from_CR(self.absZ_est[f0][V0], RC0, f0)

        # Get C from RC and R
        C = RC0 / R
    
        # Fill arrays.
        Rs[if0*len(self.V_amps[f0]) + iV0, 2] = R
        Cs[if0*len(self.V_amps[f0]) + iV0, 2] = C
    
    # To COMSOL.
    np.savetxt("COMSOL_input/Rs.txt",  np.array(Rs),   delimiter=', ', newline='\n')
    np.savetxt("COMSOL_input/Cs.txt",  np.array(Cs),   delimiter=', ', newline='\n')
  
  ############################################################################
  def create_inputs(self, opt_params=None, sim_params=None, verbose=True):
    
    self.create_RC_inputs()
    
    # Opt params to COMSOL
    if opt_params is not None:
      assert sim_params is None
      sim_params = self.params.opt_params2sim_params(opt_params)
    else:
      assert sim_params is not None
        
    for param in self.params.p_names:
      sim_params[param] *= self.get_unit(param)
      
    if verbose: print("sig={:.3g} & eps={:.3g}".format(sim_params["sigma_retina"], sim_params["epsilon_retina"]), end=', ')
        
    np.savetxt("COMSOL_input/sig.txt", np.array([[0, sim_params["sigma_retina"]],   [1, sim_params["sigma_retina"]]]),     delimiter=', ', newline='\n')
    np.savetxt("COMSOL_input/er.txt",  np.array([[0, sim_params["epsilon_retina"]], [1, sim_params["epsilon_retina"]]]),   delimiter=', ', newline='\n')

  ############################################################################
  def to_model_output(self, opt_params, rec_data):
    model_output = {}
    model_output['loss']   = self.loss.calc_loss(rec_data)
    model_output['params'] = self.params.opt_params2sim_params(opt_params)
    model_output['data']   = rec_data
    return model_output
      
  ############################################################################
  def run_parallel(self, opt_params_list, verbose=True, save_data=False, overwrite=False):
    assert len(opt_params_list) == 1
   
    if verbose: print('Clear', end=', ')
    self.clear(input=True, output=False)
    
    opt_params = opt_params_list[0]
    
    if verbose: print('Create inputs', end=', ')
    self.create_inputs(opt_params)
    
    if verbose: print('Run CMSL', end=', ')
    self.run_COMSOL(verbose=verbose)
    
    if verbose: print('Load CMSL', end=', ')
    rec_data = self.read_outputs(verbose=(verbose>=2))

    model_output_list = [self.to_model_output(opt_params, rec_data)]
    
    return model_output_list
  
  ############################################################################
  def stack_samples(self, model_output):
    if not(isinstance(model_output, list)): model_output = [model_output]
    
    optim_data = {}
    
    # Stack params.
    optim_data['params'] = {opt_param: np.concatenate([np.atleast_1d(model_output_i['params'][opt_param]) for model_output_i in model_output])\
                            for opt_param in model_output[0]['params'].keys()}
    
    # Stack loss.
    optim_data['loss']   = {loss_name: np.concatenate([np.atleast_1d(model_output_i['loss'][loss_name]) for model_output_i in model_output])\
                            for loss_name in model_output[0]['loss'].keys()}
    
    # Stack data.
    optim_data['data'] = []
    for model_output_i in model_output:
      rec_data = model_output_i['data']
        
      if isinstance(rec_data, list): optim_data['data'] += rec_data
      else:                          optim_data['data'].append(rec_data)

    return optim_data
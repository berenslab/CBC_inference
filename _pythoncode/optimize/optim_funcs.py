import pandas as pd
import numpy as np
import warnings
import time as tm
import os
import pickle
from multiprocessing import Pool
from shutil import copyfile
from matplotlib import pyplot as plt

import data_utils
   
############################################################################
############################################################################
class Optimizer():

  ############################################################################
  def __init__(
    self, params, cell=None,
    output_folder=None, t_rng=None, n_parallel=28,
    timeout=60*60*24, rec_type='optimize',
    expt_idx=0, expt_test_idx=None, test=None, check_params=True,
    raw_data_labels = ['rate BC', 'BC Vm Soma'],
    raw2model_data_labels = {'rate BC': 'rate', 'BC Vm Soma': 'Vm'},
  ):
    ''' Optimizer core.
    Calls the restim cell and manages the parameters.
    Can be combined with DELFI or similar optimization methods.
    
    Parameters:
    
    params : params object
      
    cell : cell object or None
      Cell to optimize.
      
    output_folder : str
      Folder where data will be stored in. Folder will be created in a subfolder "optim_data"
      
    t_rng : tuple
      Lower and upper bound of simulation time.
      
    n_parallel : int
      Number of parallel processes for multiprocessing.
      Should be <= the number of available CPUs.
      
    timeout : float
      Seconds after with retsim processes are terminated.
      Such simulations count as failed simulations.
      
    rec_type : str
      Recording type for retsim cell. Details are specified in retsim experiment.
      
    expt_idx : int
      Index of experiment_file in file list of cell to use as main experiment.
      
    expt_test_idx : int or None
      Index of experiment file in file list of cell to use as test experiment.
    
    test : test object or None
      Evaluates the simulation data of the test experiment.
      Must have function "is_succes"
    
    '''

    self.cell = cell
    self.params = params
    self.rec_type = rec_type
    self.n_parallel = n_parallel
    
    # Set cell? 
    if cell is not None:
      
      # Set timeout.
      self.cell.timeout = timeout
      
      # Set indexes of experiment files.
      assert len(self.cell.expt_file_list) > expt_idx
      self.expt_idx = expt_idx
      print('Main experiment: ' + self.cell.expt_file_list[expt_idx])
      
      # Set idx of test experiment file.
      if expt_test_idx is not None:
        assert len(self.cell.expt_file_list) > expt_test_idx
        assert test is not None
        print('Test experiment: ' + self.cell.expt_file_list[expt_test_idx])
        self.expt_test_idx = expt_test_idx
        self.test = test
   
      # Update t_rng.
      self.update_t_rng(t_rng=t_rng)

    # Check and set expected labels.
    self.set_expected_labels(raw_data_labels, raw2model_data_labels)
    
    # Create output folder for data storage.
    self.output_folder = output_folder
    self.init_output_folder()
    self.optim_data_file = os.path.join('optim_data', self.output_folder, 'samples.pkl')
    
    if check_params:
      self.search_params_in_cell_files(self.params, self.cell)
    
  ############################################################################
  @staticmethod
  def search_params_in_cell_files(params, cell):
    ''' Test is params are used, i.e. if they can be found in any
    retsim file of the cell.
    
    params : Parameter object
      Parameters to search.
    
    cell : Cell object
      Cell provides the file information that will be searched.
    '''

    if params is None:
      warnings.warn('Params is None. Is this correct?')
      return
      
    for param in params.p_names:
      
      if param not in params.p_range:
        print(param + ' not in range!')
      
      if param not in cell.params_default:
        print(param + ' not in default!')
          
      files_to_check = []
      files_to_check.append(os.path.join(cell.retsim_path, 'runconf', cell.chanparams_file))
      for expt_file in cell.expt_file_list:
        files_to_check.append(os.path.join(cell.retsim_path, 'expt_'+ expt_file + '.cc'))
      if cell.bp_densfile_ON  is not None:
        files_to_check.append(os.path.join(cell.retsim_path, 'runconf', cell.bp_densfile_ON))
      if cell.bp_densfile_OFF is not None:
        files_to_check.append(os.path.join(cell.retsim_path, 'runconf', cell.bp_densfile_OFF))
          
      if not data_utils.is_in_any_file(param + ';', files_to_check):
        print(param + ' not in simulation!')
 
  ############################################################################
  def set_expected_labels(self, raw_data_labels, raw2model_data_labels):
    '''Set labels that are expected from model output, and how to change them.
    Labels refer to DataFrame column names for the model output.
    
    Parameters:
    
    raw_data_labels : list of str
      Column name from model output.
      
    raw2model_data_labels : dict
      Keys (str): Raw model output names
      Values (str): Name to use instead.      
    '''
    
    assert isinstance(raw_data_labels, list)
    assert isinstance(raw2model_data_labels, dict)
    
    self.raw_data_labels = raw_data_labels
    for label in raw_data_labels:
      if label not in raw2model_data_labels.keys():
        raw2model_data_labels[label] = label
    self.raw2model_data_labels = raw2model_data_labels
 
  ############################################################################
  @staticmethod
  def copy_retsim_files(cell, retsim_folder):
    ''' Copy retsim files relevant for experiments to have a reference,
    in case the original files change.
    '''

    files_to_copy = {}

    for expt_file in cell.expt_file_list:
      file_name = 'expt_' + expt_file + '.cc'
      src_file = os.path.join(cell.retsim_path, file_name)
      trg_file = os.path.join(retsim_folder, file_name)
      files_to_copy[src_file] = trg_file
      
    for file in [cell.cone_densfile, cell.bp_densfile_ON,
                 cell.bp_densfile_OFF, cell.nval_file,
                 cell.chanparams_file]:

      if file is not None:
        src_file = os.path.join(cell.retsim_path, 'runconf', file)
        if os.path.isfile(src_file):
          trg_file = os.path.join(retsim_folder, file)
          files_to_copy[src_file] = trg_file

    for src_file, trg_file in files_to_copy.items():
      copyfile(src_file, trg_file)
      
  ############################################################################
  @staticmethod
  def check_parameter_files(cell, params, folder):
    ''' Test if cell parameters are the same as they were in a previous run.
    This is crucial if you load the data.
    '''
  
    files_vs_dicts = {}
    files_vs_dicts['cell_params_default.pkl'] = cell.params_default
    files_vs_dicts['cell_params_unit.pkl'] = cell.params_unit
    files_vs_dicts['opt_p_range.pkl'] = params.p_range

    for file, param_dict in files_vs_dicts.items():

      src_file = os.path.join(folder, file)

      if os.path.isfile(src_file):
        loaded_dict = data_utils.load_var(src_file)
        
        for param_name, param_value in param_dict.items():
          
          if param_name not in loaded_dict.keys():
            print(param_name, 'not in loaded_dict params')
          
          elif param_value != loaded_dict[param_name]:
            print(param_name, ':', param_value, '!= ', loaded_dict[param_name])

            input("Params in " + file + " are different. Press Enter to overwrite ... ")
            
      data_utils.save_var(param_dict, src_file)
              
    # If p_range was fine, p_names is fine too.
    data_utils.save_var(
      params.p_names, os.path.join(folder, 'opt_p_names.pkl')
    )
 
  ############################################################################
  def init_output_folder(self):
    ''' Create output folder.
    '''
    data_utils.make_dir('optim_data')
    assert self.output_folder is not None
    
    if os.path.exists(os.path.join('optim_data', self.output_folder)):
      print('Output folder already exists.')
    else:
      data_utils.make_dir(os.path.join('optim_data', self.output_folder))

  ############################################################################
  def get_t_rng(self):
    return self.__t_rng
  
  ############################################################################
  def get_rec_time(self):
    return self.rec_data['Time']+self.get_t_rng()[0]
  
  ############################################################################
  def update_t_rng(self, t_rng=None):
    self.__t_rng = self.cell.update_t_rng(t_rng=t_rng)
    if self.__t_rng is not None:
      self.rec_ex_size = self.estimate_rec_ex_size_from_sim_params(t_rng=self.get_t_rng(), rec_dt=self.cell.rec_dt)
    
  ############################################################################
  @staticmethod
  def estimate_rec_ex_size_from_sim_params(t_rng, rec_dt):
    return int(np.round(((t_rng[1]-t_rng[0])/rec_dt))+1)
    
  ############################################################################
  def init_rec_data(
    self, opt_params=None, sim_params=None,
    allow_loading=False, force_loading=False,
    verbose=True, save=True,
    update_ex_size=True, check_size=True,
    update_timeout=True, timout_factor=6
  ):

    # Load data?
    if allow_loading or force_loading:
      loaded_data = self.try_load_init_data(
        dirname=self.output_folder, force_loading=force_loading, verbose=verbose
      )
    else:
      loaded_data = False

    # Run experiment if data wasn't loaded.
    if not(loaded_data):
      if verbose: print('Run simulation.')
      t0 = tm.time()
      
      # Get simulation parameters.
      if (sim_params is None) and (opt_params is None):
        sim_params = self.params.p_default
      
      # Since reasonable timeout is unknown, set to high value.
      if update_timeout: self.cell.timeout = 60*60*24
      
      # Run simulation.
      self.run(
        opt_params=opt_params,
        sim_params=sim_params,
        update_cell_rec_data=True,
        update_optim_rec_data=True,
        check_ex_size=False,
        reset_retsim_stim=True,
        convert_model_output2dict=False,
        verbose=verbose,
        allow_error=True,
      );
      t_total = tm.time() - t0
      
      # Update timeout.
      if update_timeout: self.cell.timeout = timout_factor * t_total
      
      # Save init rec data.
      if save: self.save_init_data(self.output_folder)

    # Update expected recording size.
    if update_ex_size:
      self.rec_ex_size = self.rec_data['Data'].shape[0]
      print('# New expected size of recorded data: ' + str(self.rec_ex_size))
    
    # Assert recording size to match recorded data.
    # This only makes sense, when data was loaded.
    if self.rec_data['Data'].shape[0] != self.rec_ex_size:
      print('Data has size {:.2g} but should have size {:.2g}'.format(self.rec_data['Data'].shape[0], self.rec_ex_size))
    
    if verbose: self.check_rec_ex_size()
      
    if self.loss is not None:
      self.loss.update_target(rec_time=self.get_rec_time())
  
  ############################################################################
  def check_rec_ex_size(self):
    ex_size = self.estimate_rec_ex_size_from_sim_params(t_rng=self.get_t_rng(), rec_dt=self.cell.rec_dt)
    if self.rec_ex_size != ex_size:
      warnings.warn(f"Unexpected size. Does not fit t_rng and rec_dt: {ex_size} != {self.rec_ex_size}")
  
  ############################################################################
  def run(
    self, sim_params=None, opt_params=None, cell=None,
    reset_retsim_stim=True,  update_cell_rec_data=False, update_optim_rec_data=False,
    print_params=False, check_ex_size=True, verbose=False, stim_idx=None,
    convert_model_output2dict=True, allow_error=False
  ):
    ''' Run retsim with the given parameters and cell.
    '''
  
    assert not(check_ex_size and (self.rec_ex_size is None)), 'Can not check size without an expected size given'
    assert (sim_params is None) or (opt_params is None), 'Can not use both'
    
    if convert_model_output2dict: assert self.loss is not None
    
    t_start = tm.time() # Measure time.
    
    # Get retsim parameters.
    if (sim_params is None) and (opt_params is None):
      sim_params = {}
    elif opt_params is not None:
      sim_params = self.params.opt_params2sim_params(opt_params)

    # Are parameters valid.  
    params_valid = np.all([np.isfinite(p_value) for p_value in sim_params.values() if isinstance(p_value, (int, float))])
      
    # Set cell.
    cell = cell or self.cell
    assert cell is not None
      
    # Run test experiment?
    run_test = (params_valid) and (self.expt_test_idx is not None)
  
    # Run test experiment?
    if run_test:
      rec_data_test = cell.run(
        rec_type=self.rec_type, verbose=verbose, reset_retsim_stim=reset_retsim_stim,
        sim_params=sim_params, update_cell_rec_data=update_cell_rec_data,
        print_params=print_params, stim_idx=stim_idx, expt_idx=self.expt_test_idx,
      )[0]
      test_success_or_skipped = self.test.is_success(rec_data_test)
    else:
      test_success_or_skipped = True

    # Run main experiment.
    if test_success_or_skipped:
      rec_data, rec_time, rec_stim = cell.run(
        rec_type=self.rec_type, verbose=verbose, reset_retsim_stim=reset_retsim_stim,
        sim_params=sim_params, update_cell_rec_data=update_cell_rec_data,
        print_params=print_params, stim_idx=stim_idx, expt_idx=self.expt_idx
      )
    else:
      rec_data = None
      if not params_valid:
        warnings.warn('NaN parameters.')
      else:
        if verbose: warnings.warn('Test failed!')
      
    # Check rec_data.
    if allow_error:
      assert rec_data is not None, 'rec_data was None'
    
    # Check if output is as expected.
    if rec_data is not None:
      # Check data labels and rearange data to labels.
      rec_data = self.reduced_data2labels(
        self.raw_data_labels, self.raw2model_data_labels, rec_data
      )

      # Test expected size.
      if check_ex_size:
        if rec_data.shape[0] != self.rec_ex_size:
          print('Unexpected size, was {:d} but should be {:d}'.format(\
            rec_data.shape[0], self.rec_ex_size)
          )
          rec_data = None
          
      # Save data?
      if update_optim_rec_data:
        if verbose: print('Save rec_data to object.')
        self.rec_data = {'Stim': rec_stim, 'Time': rec_time, 'Data': rec_data}

    if verbose:
      if rec_data is None: print('Main experiment failed!')
      else:                print('Main experiment succeeded!')

    if convert_model_output2dict:
      return self.model_output2dict(
        sim_params=sim_params, walltime=tm.time()-t_start, rec_data=rec_data
      )
    else:
      return rec_data, rec_time, rec_stim
   
   
  ############################################################################
  @staticmethod
  def reduced_data2labels(raw_data_labels, raw2model_data_labels, rec_data):
    ''' Reduce data to labels. Take means of multiple columns.
    '''
    # Check if labels are as expected.
    if any([label not in rec_data.columns for label in raw_data_labels]):
      warnings.warn('Data labels:\n' + str(list(rec_data.columns)) +\
                    '\nExpected labels:\n' + str(raw_data_labels))
    
    # Reduce data to labels.
    if len(raw_data_labels) > 0:
      rec_data_reduced = pd.DataFrame({})
      for label in raw_data_labels:
        if rec_data[label].ndim == 1:
          rec_data_reduced[raw2model_data_labels[label]] = rec_data[label].values
        else:                         
          rec_data_reduced[raw2model_data_labels[label]] = rec_data[label].mean(axis=1).values
      
      return rec_data_reduced
    else:
      return rec_data
   
  ############################################################################
  def model_output2dict(self, sim_params, walltime, rec_data):
    ''' Out the model output into a dict with some meta-data.
    '''
    model_output = {'params': sim_params, 'wall-time': walltime}
    
    # Get labels.
    labels = None
    if len(self.raw2model_data_labels) > 0:
      labels = self.raw2model_data_labels.values()
 
    # Extract data.
    if rec_data is not None:   
      if labels is not None:
        for label in labels: model_output[label] = rec_data[label].values
      else:
        model_output['Data'] = rec_data
    
    # If no data, create NaN data.
    else:     
      if labels is not None:
        for label in labels:
          model_output[label] = np.nan*np.zeros(self.rec_ex_size)
      else:
        model_output['Data'] = None
        
    model_output['loss'] = self.loss.calc_loss(rec_data)
        
    return model_output
   
  ############################################################################
  def run_parallel(
      self, opt_params_list=None, sim_params_list=None,
      verbose=True, save_data=False, overwrite=False
    ):
    ''' Run multiple experiments in parallel.
    '''
    
    assert ((opt_params_list is None) + (sim_params_list is None)) == 1
    assert self.loss is not None
    
    if opt_params_list is not None:
      params_list = [(opt_params, True) for opt_params in opt_params_list]
    else:
      params_list = [(sim_params, False) for sim_params in sim_params_list]

    with Pool(processes=self.n_parallel) as pool:
      model_output_list = pool.starmap(self.run_parallel_instance, params_list)

    if save_data:
      self.save_samples(model_output_list, overwrite=overwrite)
    
    return model_output_list
      
  ############################################################################
  def run_parallel_instance(self, opt_or_sim_params, use_opt_params=True):
    ''' Single instance of parallel process.
    '''
    if use_opt_params:
      opt_params = opt_or_sim_params
      sim_params = None
    else:
      opt_params = None
      sim_params = opt_or_sim_params
  
    return self.run(
      opt_params=opt_params, sim_params=sim_params,
      verbose=False, reset_retsim_stim=False,
      check_ex_size=True, update_cell_rec_data=False,
      update_optim_rec_data=False, convert_model_output2dict=True
    )
  
  ############################################################################
  def save_samples(self, model_output_list, file=None, overwrite=False, verbose=False):
    ''' Save model output to file. Append data if data is already present.
    '''
    
    if file is None: file = self.optim_data_file
    
    if not overwrite:
      old_optim_data = self.load_samples(file=file, verbose=verbose)
    else:
      self.remove_samples(file=file, verbose=verbose)
      old_optim_data = None
      
    new_optim_data = self.stack_model_output_list(model_output_list)

    if old_optim_data is not None:
      optim_data = self.stack_model_output_list([old_optim_data, new_optim_data])
    else:
      optim_data = new_optim_data
    
    # Save samples.
    data_utils.save_var(optim_data, file)
  
  ############################################################################
  def remove_samples(self, file, verbose=False):   
    if os.path.exists(file):
      if verbose: print('Delete file:', file)
      os.remove(file)
  
  ############################################################################
  def load_samples(self, file=None, verbose=False):
    ''' Load samples from file.
    '''
    if file is None: file = self.optim_data_file
    
    optim_data = None
    
    if os.path.exists(file):
      if os.stat(file).st_size != 0:
        optim_data = data_utils.load_var(file)
        if verbose: print('Loaded file', file, '-> n_samples =', optim_data['loss']['total'].size)
    else:        
      if verbose: print('File', file, 'does not exist!')
    
    return optim_data
  
  ############################################################################
  def stack_model_output_list(self, model_output_list, keys_ignore_list=[]):
    ''' Stack model outputs, i.e. make arrays of parameters, loss etc.
    
    Parameters:
    
    model_output_list : list
      A list of model outputs, i.e. a list of dicts containing single samples.
      
    Returns:
    
    stacked_data : dict
      A dict containing arrays and dicts.
      Basically just a rearanged version of the input.
      Contains exactly the same data as the input.
    '''
    if not(isinstance(model_output_list, list)):
      model_output_list = [model_output_list]
    
    stacked_data = {}     

    # Concatenate params.
    stacked_data['params'] = {}
    for p_name in model_output_list[0]['params'].keys():
      stacked_data['params'][p_name] = np.concatenate([np.atleast_1d(model_output_i['params'][p_name]) for model_output_i in model_output_list])
    
    # Concatenate losses.    
    stacked_data['loss'] = {}
    for loss_name in model_output_list[0]['loss'].keys():
      stacked_data['loss'][loss_name] = np.concatenate([np.atleast_1d(model_output_i['loss'][loss_name]) for model_output_i in model_output_list])
  
    # Concatenate wall-time.
    if 'wall-time' in  model_output_list[0].keys():
      stacked_data['wall-time'] = np.concatenate([np.atleast_1d(model_output_i['wall-time']) for model_output_i in model_output_list])

    # Concatenate rest. Should be arrays.
    for key in set(model_output_list[0].keys()) - set(['params', 'loss', 'wall-time'] + keys_ignore_list):
      try:
        stacked_data[key] = np.concatenate([np.atleast_2d(model_output_i[key]) for model_output_i in model_output_list])
      except:
        warnings.warn('Could not concatenate ' + key + '. Will not be saved!')
        print("shapes to concatenate were:")
        for model_output_i in model_output_list:
          print(model_output_i[key].shape, end='\t')
        print()  

    return stacked_data
  
  ############################################################################ 
  def save_init_data(self, dirname):
    ''' Save pilot run, i.e. initial data.
    '''
    print('# Save init data to ' + dirname)
    data_utils.make_dir(os.path.join('optim_data', dirname))
    self.cell.save_init_data(os.path.join('optim_data', dirname))
    data_utils.save_var(self.rec_data,  os.path.join('optim_data', dirname, 'init_rec_data.pkl'))
      
  ############################################################################ 
  def load_init_data(self, dirname):
    ''' Load pilot run, i.e. initial data.
    '''
    self.cell.load_init_data(os.path.join('optim_data', dirname))
    self.rec_data = data_utils.load_var(os.path.join('optim_data', dirname, 'init_rec_data.pkl'))
    assert self.rec_data is not None
    
    if self.loss is not None: self.loss.update_target(rec_time=self.get_rec_time())

  ############################################################################ 
  def try_load_init_data(self, dirname, force_loading, verbose=True):
    ''' Load data if possible and return True if succeded.
    '''
    # Load?
    if force_loading:
      self.load_init_data(dirname)
      success = True
    else:
      try:
        self.load_init_data(dirname)
        success = True
      except:
        success = False
        
    if verbose:
      if success: print('Loaded successfully.')
      else:       print('Loading failed!')

    return success
    
############################################################################
############################################################################
class EmptyOptimizer(Optimizer):

  def __init__(self, *args, **kwargs):
    pass
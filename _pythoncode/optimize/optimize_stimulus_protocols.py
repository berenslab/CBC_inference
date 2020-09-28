import numpy as np
from pandas import DataFrame
from matplotlib import pyplot as plt
from time import sleep
from multiprocessing import Pool
import os
import time

import comsol_utils
import data_utils
from optim_funcs import Optimizer

DBP1 = 6
HBP1 = 10

global comsol2retsim_folder
comsol2retsim_folder = os.path.join(
  '/gpfs01', 'berens', 'user', 'joesterle', 'berens', 'project_bipolarcells', 'COMSOL2retsim_interface/'
)

############################################################################
def write_to_file(data, file):
  with open(file, 'w+') as f:
    f.write(str(data))
  assert os.path.isfile(file)

############################################################################
class OptimizerStimulus(Optimizer): 

  def __init__(
      self, cells, output_folder, t_rng,
      comsol_filename, comsol_batch_size,
      comsol_samples_prefix, comsol_samples_suffix,
      params, comsol_global_parameters={},
      verbose=True, rec_type='optimize', 
      n_reps=1, inter_rep_dt=0.0, predur_stim=0.0,
      set_comsol2retsim_folder=None,
      **kwargs
    ):
    
    self.DEBUG = False
    global comsol2retsim_folder
    
    if set_comsol2retsim_folder:
      comsol2retsim_folder = set_comsol2retsim_folder
    
    # Initialize optimizer.
    super().__init__(
      cell=None, params=None, rec_type=rec_type, check_params=False,
      output_folder=output_folder, **kwargs
    )
   
    self.params = params
    
    # Set cells and cell target.
    self.cells = cells
    for cell in self.cells: cell.stim = None
    
    # Set sim dt.
    self.sim_dt = self.cells[0].sim_dt
    for cell in self.cells[1:]: assert self.sim_dt == cell.sim_dt
      
    # Set stim time.
    self.__t_rng = t_rng
    self.stim_time = np.arange(self.__t_rng[0],self.__t_rng[1],self.sim_dt)
    
    # Set stimulus repetition parameters.
    self.n_reps = n_reps
    self.inter_rep_dt = inter_rep_dt
    self.predur_stim = predur_stim
    
    for cell in cells: cell.params_default['rec_predur'] = predur_stim
    
    # Update COMSOL parameters.
    self.set_comsol_filename(comsol_filename)
    self.set_samples_prefix_and_suffix(comsol_samples_prefix, comsol_samples_suffix)
    self.set_comsol_batch_size(comsol_batch_size)
    self.set_global_parameters(comsol_global_parameters)
    
    # Delete all stimuli in folder.
    self.set_input_not_ready()
    self.wait_for_comsol_idle()
    self.set_output_not_ready()
    self.clean_input_samples()
  
  ############################################################################
  @staticmethod
  def set_comsol_batch_size(batch_size):
    file = os.path.join(comsol2retsim_folder , 'interface', '_batch_size')
    write_to_file(batch_size, file)
  
  ############################################################################
  @staticmethod
  def set_comsol_filename(comsol_filename):
    file = os.path.join(comsol2retsim_folder, 'interface', '_comsol_filename')
    write_to_file(comsol_filename, file)
      
  ############################################################################
  def set_samples_prefix_and_suffix(self, prefix, suffix):
    prefixfile = os.path.join(comsol2retsim_folder, 'interface', '_samples_prefix')
    write_to_file(prefix, prefixfile)
    self.comsol_samples_prefix = prefix
      
    suffixfile = os.path.join(comsol2retsim_folder, 'interface', '_samples_suffix')
    write_to_file(suffix, suffixfile)
    self.comsol_samples_suffix = suffix
  
  ############################################################################
  @staticmethod
  def set_global_parameters(comsol_global_parameters):
  
    for k, v in comsol_global_parameters.items():
      outputfile = os.path.join(comsol2retsim_folder, 'comsol_input', 'global', k  + ".txt")
      
      if isinstance(v, (float, int)):
        np.savetxt(outputfile, np.array([[0, v], [1, v]]), delimiter=', ', newline='\n')
      elif isinstance(v, np.ndarray, list):
        np.savetxt(outputfile, np.asarray(v), delimiter=', ', newline='\n')
      else:
       raise TypeError
  
  ############################################################################
  @staticmethod
  def clean_input_samples():
    files = os.listdir(os.path.join(comsol2retsim_folder, 'comsol_input', 'samples'))
    for file in files:
      os.remove(os.path.join(comsol2retsim_folder, 'comsol_input', 'samples', file))
  
  ############################################################################
  @staticmethod
  def wait_for_comsol_idle( verbose=True):
    if verbose: print('Waiting for COMSOL idle ... ', end='')
    
    file = os.path.join(comsol2retsim_folder, 'interface', '_comsol_idle')
    
    while (True):
      sleep(0.1)
      
      with open(file, 'r') as f:
        line = f.readline()
      
      if line == '1':
        break
    
    if verbose: print('Ready!')

  ############################################################################
  @staticmethod
  def set_input_ready():
    file = os.path.join(comsol2retsim_folder, 'interface', '_input_ready')
    write_to_file(1, file)

  ############################################################################
  @staticmethod      
  def set_input_not_ready():
    file = os.path.join(comsol2retsim_folder, 'interface', '_input_ready')
    if os.path.isfile(file): os.remove(file)
      
  ############################################################################
  @staticmethod      
  def set_output_not_ready():
    file = os.path.join(comsol2retsim_folder, 'interface', '_output_ready')
    if os.path.isfile(file): os.remove(file)
  
  ############################################################################
  def run_comsol(self, verbose=True):
      
    # Set input ready.
    self.set_input_ready()
    
    # Wait for output being ready.
    output_ready = False
    if verbose: print('Waiting for COMSOL ... ', end='')
    
    while (True):
      if '_output_ready' in os.listdir(os.path.join(comsol2retsim_folder, 'interface')):
        self.set_output_not_ready()
        break
      
      sleep(0.01)
    
    if verbose: print('Ready!')
      
  ############################################################################
  def set_stim_generator(self, stim_generator):
    self.stim_generator = stim_generator
    
  ############################################################################
  def set_loss(self, loss):
    self.loss = loss

  ############################################################################
  def set_comsol_input(self, params_list):
    
    assert self.stim_generator is not None
    
    self.clean_input_samples()
    
    # Create stimuli for all parameters.
    for idx, params in enumerate(params_list):
      stimulus_filename = os.path.join(
        comsol2retsim_folder, 'comsol_input', 'samples',
        self.comsol_samples_prefix + str(idx) + self.comsol_samples_suffix
      )
      self.stim_generator.create_stimulus(params=params, filename=stimulus_filename)
      
  ############################################################################
  def set_init_comsol_input(self):
    
    assert self.stim_generator is not None
    
    self.clean_input_samples()
    
    # Create stimuli for all parameters.
    stimulus_filename = os.path.join(
      comsol2retsim_folder, 'comsol_input', 'samples',
      self.comsol_samples_prefix + str(0) + self.comsol_samples_suffix
    )
    self.stim_generator.create_zero_stimulus(filename=stimulus_filename)
  
  ############################################################################
  def get_comsol_output(self, n_stimuli):
  
    # Set output to not be ready anymore.
    self.set_output_not_ready()
  
    for cell in self.cells: cell.stim = None
  
    Vext_files = sorted(os.listdir(comsol2retsim_folder + "comsol_output/samples/"))
    
    for cell in self.cells:
      if cell.bp_type_idx == DBP1:
        Vext_files_cell = ['Vext_ON_BC_idx_'  + str(stim_idx) + '.csv' for stim_idx in range(n_stimuli)]
      elif cell.bp_type_idx == HBP1:
        Vext_files_cell = ['Vext_OFF_BC_idx_' + str(stim_idx) + '.csv' for stim_idx in range(n_stimuli)]
      else:
        raise
      
      # Check that all files are present.
      for Vext_file in Vext_files_cell:
        assert Vext_file in Vext_files, (Vext_file + ' not in ' + str(Vext_files))
      
      # Set stimuli
      cell_stims = []
      for Vext_file in Vext_files_cell:
        stim_file = comsol2retsim_folder + "comsol_output/samples/" + Vext_file
        loaded_stim = comsol_utils.comsol2dataframe(stim_file)
        cell_stims.append(loaded_stim)
        
      cell.set_stim(cell_stims, n_reps=self.n_reps, inter_rep_dt=self.inter_rep_dt)
  
  ############################################################################
  def run_parallel(
      self, params_list, verbose=False, save_data=False, overwrite=False, 
      skip_comsol=False,
    ):
      
    if not skip_comsol:
      self.set_comsol_input(params_list)
      self.run_comsol(verbose=verbose)
      self.get_comsol_output(n_stimuli=len(params_list))
    
    # Run all cells and stimuli in parallel if possible.
    params_list_parallel = []
    for stim_idx in range(len(params_list)):
      for cell in self.cells:
        params_list_parallel.append([stim_idx, cell])
        
    with Pool(processes=self.n_parallel) as pool:
      rec_data_list = pool.starmap(self.run_instance, params_list_parallel)

    if self.DEBUG: return rec_data_list

    # Rearrange data.
    n_per_stim = len(self.cells)
    rec_data_list = [rec_data_list[stim_idx*n_per_stim:(stim_idx+1)*n_per_stim] for stim_idx in range(len(params_list))]
    
    assert len(rec_data_list) == len(params_list)
   
    model_output_list = []
    for params, rec_data in zip(params_list, rec_data_list):
      assert len(rec_data) == n_per_stim
      model_output = {cell.bp_type: {rec_type: rec_data_i['Data'][rec_type].values for rec_type in ['rate', 'Vm']} for rec_data_i, cell in zip(rec_data, self.cells)}
      model_output.update({'params': params, 'loss': self.loss.calc_loss(rec_data_dict=model_output)})
      model_output_list.append(model_output)
    
    if save_data: self.save_samples(model_output_list, overwrite=overwrite)
    
    return model_output_list 

  ############################################################################
  def run_instance(self, stim_idx, cell, overwrite=False, verbose=False):
    
    # Simulate cell.
    assert self.rec_ex_size is not None, 'Initialize first'
    assert cell is not None
    
    # Run simulation.
    rec_data, rec_time, rec_stim = self.run(
      cell=cell, stim_idx=stim_idx, sim_params={}, verbose=verbose,
      reset_retsim_stim=True, check_ex_size=True,
      update_cell_rec_data=False, update_optim_rec_data=False,
      convert_model_output2dict=False
    )

    # Discard all but the last rep.
    rec_data = rec_data.iloc[int(rec_data.shape[0]*((self.n_reps-1)/self.n_reps)):,:]

    # Discard predur.
    idx = np.argmax(np.array(rec_time) >= self.predur_stim)
    rec_data = rec_data.iloc[idx:,:]
    
    return {'stim_idx': stim_idx, 'cell': cell.bp_type, 'Data': rec_data}

  ############################################################################
  def stack_model_output_list(self, model_output_list, keys_ignore_list=None):
    if not(isinstance(model_output_list, list)):
      model_output_list = [model_output_list]
    
    optim_data = {}
    
    # Concatenate params.
    optim_data['params'] = np.vstack([np.atleast_1d(model_output['params']) for model_output in model_output_list])
    
    # Concatenate losses.    
    optim_data['loss'] = {}
    for loss_name in model_output_list[0]['loss'].keys():
      optim_data['loss'][loss_name] = np.concatenate([np.atleast_1d(model_output['loss'][loss_name]) for model_output in model_output_list])
    
    # Concatenate cell data.
    for cell in self.cells:
      optim_data[cell.bp_type] = {}
      optim_data[cell.bp_type]['rate'] = np.vstack([np.atleast_1d(model_output[cell.bp_type]['rate']) for model_output in model_output_list])
      optim_data[cell.bp_type]['Vm']   = np.vstack([np.atleast_1d(model_output[cell.bp_type]['Vm'])   for model_output in model_output_list])

    return optim_data
      
  ############################################################################
  def init_rec_data(
      self, allow_loading=False, force_loading=False, verbose=True,
      **kwargs
    ):
    
    for cell in self.cells:
      assert not cell.make_cones
    
    # Load data?
    if allow_loading or force_loading:
      loaded_data = self.try_load_init_data(
        dirname=self.output_folder, force_loading=force_loading, verbose=verbose
      )
    else:
      loaded_data = False
      
    # Initialize?
    if not(loaded_data):
      self.rec_data = self.run_init_rec_data(verbose=verbose, **kwargs)
      self.save_init_data(self.output_folder)

    # Set rec_size.
    self.set_rec_ex_size(verbose=verbose)
    
    # Save stimulation time.
    data_utils.save_var(self.stim_time, file='optim_data/' + self.output_folder + '/stim_time.pkl')
  
  ############################################################################ 
  def run_init_rec_data(self, verbose, eq_exists=False):
  
    # Initialize with zero stimulus.
    self.set_init_comsol_input()
    self.run_comsol()
    self.get_comsol_output(n_stimuli=1)
    
    if verbose: print('Compute cell responses:')
    
    # Prepare cells.
    for cell in self.cells:
      cell.params_default.update({"save_eq": 1, "load_eq": 0, "run_predur_only": 1})
        
    # Run to equilibrium.
    if not eq_exists:
      if verbose: print("Create equilibrium files.")
      with Pool(processes=self.n_parallel) as pool:
        pool.map(self.run_init_rec_data_instance_eq, self.cells);

    # Set to loading.
    for cell in self.cells:
      cell.params_default.update({"save_eq": 0, "load_eq": 1, "run_predur_only": 0})   

    # Run simulation in final state.
    if verbose: print("Run experiments after loading equilibrium files.")
    with Pool(processes=self.n_parallel) as pool:
      rec_data_list = pool.map(self.run_init_rec_data_instance_load, self.cells)
        
    # Save combined for all cells.
    rec_data = {cell.bp_type: {'Data': rec_data[0], 'Time': rec_data[1]} for cell, rec_data in zip(self.cells, rec_data_list)}
    return rec_data
  
  ###########################################################################
  def run_init_rec_data_instance_eq(self, cell):
    return self.run(
      cell=cell, update_cell_rec_data=False, update_optim_rec_data=False,
      check_ex_size=False, stim_idx=0, convert_model_output2dict=False,
      verbose=True, print_params=False,
    )
  
  ###########################################################################
  def run_init_rec_data_instance_load(self, cell):
    return self.run(
      cell=cell, update_cell_rec_data=True, update_optim_rec_data=True,
      check_ex_size=False, stim_idx=0, convert_model_output2dict=False,
      verbose=True, print_params=False,
    )
  
  ###########################################################################
  def set_rec_ex_size(self, verbose=True):
    assert self.rec_data is not None
    rec_ex_size_list = [self.rec_data[cell.bp_type]['Time'].size for cell in self.cells]
    if not np.unique(rec_ex_size_list).size == 1:
      print('Warning: Not all cell outputs should have the same size.')
    self.rec_ex_size = rec_ex_size_list[0]    
    if verbose: print('Expected output size of simulation: ' + str(self.rec_ex_size))
  
  ############################################################################ 
  def save_init_data(self, dirname):
    data_utils.make_dir('optim_data/' + dirname)
    for cell in self.cells:
      data_utils.make_dir('optim_data/' + dirname + '/' + cell.bp_type)
      cell.save_init_data('optim_data/' + dirname + '/' + cell.bp_type)
    data_utils.save_var(self.rec_data, file='optim_data/' + dirname + '/init_rec_data.pkl')

  ############################################################################ 
  def load_init_data(self, dirname):
    for cell in self.cells:
      cell.load_init_data('optim_data/' + dirname + '/' + cell.bp_type)
    self.rec_data = data_utils.load_var('optim_data/' + dirname + '/init_rec_data.pkl')
    assert self.rec_data is not None    
    
  ############################################################################    
  def plot_samples(self, samples, n_max=100):
    if samples is None: print('No samples given!'); return
    n_samples = samples['loss']['total'].size
    # Plot.
    plt.figure(figsize=(12,5))
    plt.subplot(311)
    plt.ylabel('Vesicles/s')
    plt.xlabel('Time [s]')
    if n_samples > n_max: print('Plot ' + str(n_max) + ' random samples!')
    idx = np.random.choice(np.arange(0,n_samples), size=np.min([n_max, n_samples]), replace=False)
    for cell_idx, cell in enumerate(self.cells):
      plt.plot(self.rec_data[cell.bp_type]['Time'], samples[cell.bp_type]['rate'][idx,:].T, c=['r', 'b'][cell_idx])
    ax = plt.subplot(312)
    plt.ylabel('Stim [uA]')
    for idx_i in idx:
      current = self.stim_generator.create_stimulus(samples['params'][idx_i,:], plot=False, filename=None)
      plt.plot(self.stim_time, current)
    plt.ylim((-np.max(np.abs(current)), np.max(np.abs(current))))
    plt.subplot(313)
    plt.ylabel('Loss')
    plt.plot(samples['loss']['total'], '.-', label='total', linewidth=0.3, markersize=1)
    plt.xlabel('# samples')
    plt.ylabel('Loss')
    plt.legend(bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()
    
  ############################################################################    
  def plot_init_rec_data(self):
  
    # Plot trace.
    fig, axs = plt.subplots(2,1,figsize=(12,4), sharex=True)
    for cell in self.cells:
      axs[0].plot(self.rec_data[cell.bp_type]['Time'], self.rec_data[cell.bp_type]['Data']['Vm'], label=cell.bp_type)
      axs[1].plot(self.rec_data[cell.bp_type]['Time'], self.rec_data[cell.bp_type]['Data']['rate'], label=cell.bp_type)
    axs[0].legend(loc='upper right')
    axs[0].set_ylabel('Vm')
    axs[1].set_ylabel('rate')
    axs[1].set_xlabel('Time [ms]')
    
    plt.tight_layout()
    plt.show()
    
############################################################################
class OptimizerStimulusMultiParams(OptimizerStimulus): 

  ############################################################################
  def __init__(
      self, *args, cell_params_list, **kwargs,
    ):
    
    # Initialize optimizer.
    super().__init__(*args, **kwargs)
    self.cell_params_list = cell_params_list
    
    # Update cell params to load different eq files.
    for cell, cell_params in zip(self.cells, self.cell_params_list):
      for i, cell_params_i in enumerate(cell_params):
        if 'eqfile' not in cell_params_i.keys():
          cell_params_i['eqfile'] = cell.bp_type + '_params_idx_'+ str(i) + '.eq' 
    
  ############################################################################
  def run_parallel(
      self, params_list, verbose=False, save_data=False, overwrite=False,
      skip_comsol=False
    ):
      
    assert self.loss is not None
      
    if not skip_comsol:
      self.set_comsol_input(params_list)
      self.run_comsol(verbose=verbose)
      self.get_comsol_output(n_stimuli=len(params_list))
    
    # Run all cells and stimuli in parallel if possible.
    params_list_for_parallel_computation = []
    for stim_idx in range(len(params_list)):
      for cell, cell_params in zip(self.cells, self.cell_params_list):
        cell.create_retsim_stim_file(stim_idx=stim_idx)
        for cell_params_i in cell_params:
          params_list_for_parallel_computation.append([stim_idx, cell, cell_params_i])
        
    with Pool(processes=self.n_parallel) as pool:
      rec_data_list = pool.starmap(self.run_instance, params_list_for_parallel_computation)
    
    # Rearrange data.
    n_per_stim = int(len(rec_data_list) / len(params_list))
    rec_data_list = [rec_data_list[stim_idx*n_per_stim:(stim_idx+1)*n_per_stim] for stim_idx in range(len(params_list))]

    if self.DEBUG: return rec_data_list
    
    model_output_list = []
    for params, rec_data in zip(params_list, rec_data_list):
      assert len(rec_data) == n_per_stim
      model_output = {}
      for i, (cell, cell_params) in enumerate(zip(self.cells, self.cell_params_list)):
        model_output[cell.bp_type] = {rec_type: np.full((1,self.rec_ex_size,len(cell_params)), np.nan) for rec_type in ['rate', 'Vm']}
        
        rec_data_cell = rec_data[i*len(cell_params):(i+1)*len(cell_params)] 
        for j, rec_data_cell_i in enumerate(rec_data_cell):
          for rec_type in ['rate', 'Vm']:
            model_output[cell.bp_type][rec_type][0,:,j] = rec_data_cell_i['Data'][rec_type].values
      
      model_output.update({'params': params, 'loss': self.loss.calc_loss(rec_data_dict=model_output)})
      model_output_list.append(model_output)
      
    assert len(model_output_list) == len(params_list)
    
    if save_data: self.save_samples(model_output_list, overwrite=overwrite)
    
    return model_output_list 

  ############################################################################
  def run_instance(self, stim_idx, cell, cell_params, overwrite=False, verbose=False):
    
    # Simulate cell.
    assert self.rec_ex_size is not None, 'Initialize first'
    assert cell is not None
    
    # Run simulation.
    rec_data, rec_time, rec_stim = self.run(
      cell=cell, stim_idx=stim_idx, sim_params=cell_params, verbose=verbose,
      reset_retsim_stim=False, check_ex_size=True,
      update_cell_rec_data=False, update_optim_rec_data=False,
      convert_model_output2dict=False
    )

    # Discard all but the last rep.
    rec_data = rec_data.iloc[int(rec_data.shape[0]*((self.n_reps-1)/self.n_reps)):,:]

    # Discard predur.
    idx = np.argmax(np.array(rec_time) >= self.predur_stim)
    rec_data = rec_data.iloc[idx:,:]
    
    return {'stim_idx': stim_idx, 'cell': cell.bp_type, 'Data': rec_data, 'cell_params': cell_params}
    
  ############################################################################ 
  def run_init_rec_data(self, verbose, eq_exists=False):
  
    # Initialize with zero stimulus.
    self.set_init_comsol_input()
    self.run_comsol()
    self.get_comsol_output(n_stimuli=1)
    
    if verbose: print('Compute cell responses:')
    
    parallel_params_list = []
    for cell, cell_params in zip(self.cells, self.cell_params_list):
      for cell_params_i in cell_params:
        parallel_params_list.append((cell, cell_params_i))
    
    # Prepare cells.
    if not eq_exists:
      for cell in self.cells:
        cell.params_default.update({"save_eq": 1, "load_eq": 0, "run_predur_only": 1})
          
      # Run to equilibrium.
      if verbose: print("Create equilibrium files.")
      with Pool(processes=self.n_parallel) as pool:
        pool.starmap(self.run_init_rec_data_instance_eq, parallel_params_list);

    # Set to loading.
    for cell in self.cells:
      cell.params_default.update({"save_eq": 0, "load_eq": 1, "run_predur_only": 0})

    # Run simulation in final state.
    if verbose: print("Run experiments after loading equilibrium files.")
    with Pool(processes=self.n_parallel) as pool:
      rec_data_list = pool.starmap(self.run_init_rec_data_instance_load, parallel_params_list)
        
    # Save combined for all cells.
    for rec_data_i in rec_data_list:
      assert np.all(rec_data_i[1] == rec_data_list[0][1])
      assert rec_data_i[0].shape == rec_data_list[0][0].shape
    
    rec_data = {cell.bp_type: {'Time': rec_data_list[0][1], 'Data': []} for cell in self.cells}
    
    for (cell, *_), rec_data_i in zip(parallel_params_list, rec_data_list):
      rec_data[cell.bp_type]['Data'].append(rec_data_i[0])
    
    return rec_data
  
  ###########################################################################
  def run_init_rec_data_instance_eq(self, cell, cell_params):
    print('Create eq file:', cell.bp_type, '-->', cell_params['eqfile'])
    return self.run(
      cell=cell, sim_params=cell_params,
      update_cell_rec_data=False, update_optim_rec_data=False,
      check_ex_size=False, stim_idx=0, convert_model_output2dict=False,
      verbose=False, print_params=False,
    );
  
  ###########################################################################
  def run_init_rec_data_instance_load(self, cell, cell_params):
    print('Run with eq file:', cell.bp_type, '-->', cell_params['eqfile'])
    return self.run(
      cell=cell,  sim_params=cell_params,
      update_cell_rec_data=True, update_optim_rec_data=True,
      check_ex_size=False, stim_idx=0, convert_model_output2dict=False,
      verbose=False, print_params=False,
    );
    
  ############################################################################    
  def plot_init_rec_data(self):
    fig, axs = plt.subplots(2,1,figsize=(12,4), sharex=True)
    for cell, color in zip(self.cells, ['r', 'b']):
      for i, rec_data_i in enumerate(self.rec_data[cell.bp_type]['Data']):
        axs[0].plot(self.rec_data[cell.bp_type]['Time'], rec_data_i['Vm'], label=cell.bp_type if i == 0 else '_', c=color)
        axs[1].plot(self.rec_data[cell.bp_type]['Time'], rec_data_i['rate'], label=cell.bp_type if i == 0 else '_', c=color)
    axs[0].legend(loc='upper right')
    axs[0].set_ylabel('Vm')
    axs[1].set_ylabel('rate')
    axs[1].set_xlabel('Time [ms]')
    
    plt.tight_layout()
    plt.show()
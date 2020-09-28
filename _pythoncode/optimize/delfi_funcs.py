import numpy as np
import time as tm
import os
import warnings
import importlib
from shutil import copyfile
from copy import deepcopy

import optim_funcs
import data_utils
import my_delfi_funcs

############################################################################
############################################################################
class DELFI_Optimizer():
  
  nn_epochs = 100
  nn_minibatch = 100
  
  ############################################################################
  def __init__(
    self, optim, prior, n_parallel=25, gen_minibatch=300,
    samples_file=None, snpe_type='B', round_cl=1, scalar_loss=True,
    post_as_truncated_normal=True,
    samples_folder='samples', snpe_folder='snpe',
    backups_folder='backups', retsim_folder='retsim',
  ):
    ''' Initialize DELFI wrapper.
    
    Parameters:
    
    post_as_truncated_normal : bool
      Transform posteriors to truncated normal distributions.
    '''
  
    self.snpe_type = snpe_type
    self.prior = prior
   
    self.optim = optim
    self.optim.n_parallel = n_parallel
    self.scalar_loss = scalar_loss
      
    self.gen_minibatch = gen_minibatch
    
    self.round_cl = round_cl
    
    # Make folders to store data.
    self.general_folder = os.path.join('optim_data', self.optim.output_folder)
    self.samples_folder = os.path.join(self.general_folder, samples_folder+'/')
    self.snpe_folder = os.path.join(self.general_folder, snpe_folder+'/')
    self.backups_folder = os.path.join(self.general_folder, backups_folder+'/')
    self.retsim_folder = os.path.join(self.general_folder, retsim_folder+'/')
    
    data_utils.make_dir(self.samples_folder)
    data_utils.make_dir(self.snpe_folder)
    data_utils.make_dir(self.backups_folder)
    data_utils.make_dir(self.retsim_folder)
  
    # Save retsim files and check parameters.
    if self.optim.cell is not None:
      self.optim.copy_retsim_files(
        cell=self.optim.cell, retsim_folder=self.retsim_folder
      )
      self.optim.check_parameter_files(
        cell=self.optim.cell, params=self.optim.params, folder=self.general_folder
      )

    # Define file for samples.
    self.update_samples_file(samples_file)
  
    self.post_as_truncated_normal = post_as_truncated_normal
    
    if self.post_as_truncated_normal:
      self.trunc_lower = self.prior.lower
      self.trunc_upper = self.prior.upper
    else:
      self.trunc_lower = None
      self.trunc_upper = None
  
  ############################################################################
  def create_delfi_generator(self):
    '''Initialize DELFI functions.
    
    Parameters:
    
    lower, upper : 1d-arrays with size of parameters
      Lower and upper bounds of parameters.
      Used in generator.
      Parameters outside this range will be resampled.
    '''
  
    delfi_m = my_delfi_funcs.retsim_Simulator(
      self, dim_param=self.optim.params.p_N
    )
    
    if self.scalar_loss:
      n_summary = 1
    else:
      n_summary = len(self.optim.loss.loss_params)
      
    delfi_s = my_delfi_funcs.retsim_SummaryStats(
      n_summary=n_summary,
    )
    
    delfi_generator = my_delfi_funcs.retsim_Generator(
      model=delfi_m, prior=self.prior, summary=delfi_s,
      gen_minibatch=self.gen_minibatch
    )
    
    return delfi_generator
  
  ############################################################################
  def init_SNPE(
    self, verbose=False,
    n_components=1,
    obs=None,
    svi=False,
    pseudo_obs_dim=0,
    pseudo_obs_n=None,
    pseudo_obs_perc=None, 
    kernel_bandwidth=2,
    kernel_bandwidth_n=None,
    kernel_bandwidth_perc=None,
    kernel_bandwidth_min=0.0,
    pseudo_obs_use_all_data=False,
    prior_mixin=0.0,
    use_doubling=False,
    seed=777,
    **kwargs,
  ):
  
    ''' Initialize SNPE solver.
    '''
  
    from delfi import inference
  
    np.random.seed(seed)
    self.prior.reseed(seed)    
    self.save_random_state()

    from delfi import __version__ as delfi_version
    n_params = self.optim.params.p_N
    
    print(f'Create new network, using SNPE-{self.snpe_type}')
    
    if obs is None:
      n_loss = 1 if self.scalar_loss else len(self.optim.loss.loss_params)
      obs = np.zeros((1,n_loss))

    # SNPE-B, modified.
    if self.snpe_type in ['b', 'B']:

      assert delfi_version == '0.5.1', 'Use the self implemented version if you want to use this function.'
      
      if isinstance(kernel_bandwidth, float):
        kernel_bandwidth = np.full(obs.size, kernel_bandwidth)
      
      from StrechedGauss import StrechedGauss
      kern = StrechedGauss(obs=obs, bandwidth=kernel_bandwidth, max_weight_range=0.0)

      n_hiddens = [int(((n_params*n_params)/2 + n_params))*n_components]*3
      print("n_hiddens:", n_hiddens)

      self.inf_snpe = inference.SNPE(
        generator=self.create_delfi_generator(),
        n_components=n_components,
        obs=obs,
        pseudo_obs_dim=pseudo_obs_dim,
        pseudo_obs_n=pseudo_obs_n,
        pseudo_obs_perc=pseudo_obs_perc,
        kernel_bandwidth_n=kernel_bandwidth_n,
        kernel_bandwidth_perc=kernel_bandwidth_perc,
        kernel_bandwidth_min=kernel_bandwidth_min,
        pseudo_obs_use_all_data=pseudo_obs_use_all_data,
        n_hiddens=n_hiddens,
        pilot_samples=None,
        convert_to_T=None,
        prior_norm=False,
        prior_mixin=prior_mixin,
        verbose=verbose,
        svi=svi,
        kernel=kern,
        n_inputs=(obs.size,),
        n_outputs=self.optim.params.p_N,
        use_doubling=use_doubling,
      )
      
    # SNPE-C, modified.
    if self.snpe_type in ['c', 'C']:

      assert delfi_version != '0.5.1', 'Use the latest version if you want to use this function.'

      self.inf_snpe = inference.APT(
        generator=self.create_delfi_generator(),
        obs=obs,
        n_hiddens=[int(((n_params*n_params)/2 + n_params) * n_components)]*3,
        pilot_samples=0,
        prior_norm=False,
        verbose=verbose,
        **kwargs,
      )
    
  ############################################################################
  def save_random_state(self):
    ''' Save numpy random state to folder.
    '''
    data_utils.save_var(np.random.get_state(), os.path.join(self.general_folder, 'random_state.pkl'))
    
  ############################################################################
  def load_random_state(self):
    ''' Load numpy random state from folder.
    '''
    np.random.set_state(data_utils.load_var(os.path.join(self.general_folder, 'random_state.pkl')))
    
  ############################################################################
  def run_SNPE(
      self, n_samples_per_round,
      max_duration_minutes=np.inf, max_rounds=None,
      continue_optimization=False, load_init_tds=False,
    ):
    
    ''' Run SNPE optimization.
    Add some utility functionality to the vanilla SNPE:
    Store all relevant data to folder.
    Optimization can be continued later, loading those files.
    One can also use only the samples from the prior and rerun the rest.
    
    Parameters:
    
    n_samples_per_round : int
      Number of samples generated per NN training round.
      
    max_duration_minutes : float
      If this time is exceeded, will stop after the next round.
      
    max_rounds : int
      Maximum number of rounds.
      
    continue_optimization : bool
      Load previously produced data and continue.
      
    load_init_tds : bool
      Load initial training data, i.e. will use old training data of last round.
      Otherwise algorithm will start from strech.
      Can be used to play around with NN parameters.
 
    '''
    
    assert not(load_init_tds and continue_optimization)

    self.load_random_state()

    #### Initialize.
    if continue_optimization:
      inf_snpes, logs, tds, sample_distributions, n_samples, \
        kernel_bandwidths, pseudo_obs = self.load_SNPE_rounds()
      self.inf_snpe = inf_snpes[-1]
    else:
      self.create_SNPE_rounds_backup()
    
      assert self.inf_snpe is not None, 'Initialize network first'
      inf_snpes            = [deepcopy(self.inf_snpe)]
      logs                 = []
      tds                  = []
      sample_distributions = [self.prior]
      n_samples            = []     
      if self.snpe_type in ['b', 'B']:
        kernel_bandwidths  = []
        pseudo_obs         = []
    
    ### Load data?
    load_tds = False    
    if load_init_tds:
      load_tds = True
      # Restore file from backups.
      copyfile(os.path.join(self.backups_folder, 'delfi_samples_r0.pkl'),
               os.path.join(self.samples_folder + 'delfi_samples_r0.pkl'))
      load_file = os.path.join(self.samples_folder, 'delfi_samples_r0.pkl')

    elif continue_optimization:
      files = sorted(os.listdir(self.samples_folder))
      if len(n_samples) < len(files):
        print('Found incomplete round. Will load samples.')
        load_tds = True
        load_file = self.samples_folder + files[-1]
    
    if load_tds:
      loaded_tds, n_loaded_samples = self.load_tds_from_file(
        file=load_file, params=self.optim.params
      )
    else:
      loaded_tds = None
      n_loaded_samples = 0
      
    ### Optimize.
    if max_rounds is None: max_rounds = np.inf
    t0 = tm.time()
    i_round = 0
    
    print()
    
    while (((tm.time() - t0)/60) < max_duration_minutes) and (self.inf_snpe.round < max_rounds):
      if i_round > 0:
        loaded_tds = None
        n_loaded_samples = 0
              
      # Update output file.
      samples_file = self.samples_folder + 'delfi_samples_r{:d}.pkl'.format(self.inf_snpe.round)
      self.update_samples_file(samples_file)
      
      print('# Sample output file:', self.samples_file)
      print('# Round', self.inf_snpe.round+1, '- Running for {:.2f} [min]'.format(((tm.time()-t0)/60)))
      
      ### Run inference method.     
      new_log, new_td, new_posterior = self.run_SNPE_round(
        n_train=[n_samples_per_round], proposal=sample_distributions[-1],
        initial_tds=loaded_tds, verbose=True
      )
      
      if self.post_as_truncated_normal:
        new_posterior = my_delfi_funcs.normal2truncated_normal(
          new_posterior, lower=self.trunc_lower, upper=self.trunc_upper
        ) 
      
      ### Append data.
      logs.append(new_log)
      sample_distributions.append(new_posterior)
      tds.append(new_td)
      inf_snpes.append(deepcopy(self.inf_snpe))
      n_samples.append(new_td[1].shape[0])
      if self.snpe_type in ['b', 'B']:
        pseudo_obs.append(self.inf_snpe.pseudo_obs[-1].copy())
        kernel_bandwidths.append(self.inf_snpe.kernel_bandwidth[-1].copy())
      
      ### Save data to files.
      self.save_random_state()
      data_utils.save_var(inf_snpes,            file=os.path.join(self.snpe_folder, 'inf_snpes.pkl'))
      data_utils.save_var(sample_distributions, file=os.path.join(self.snpe_folder, 'sample_distributions.pkl'))
      data_utils.save_var(logs,                 file=os.path.join(self.snpe_folder, 'logs.pkl'))
      data_utils.save_var(tds,                  file=os.path.join(self.snpe_folder, 'tds.pkl'))
      data_utils.save_var(n_samples,            file=os.path.join(self.snpe_folder, 'n_samples.pkl'))
      if self.snpe_type in ['b', 'B']:
        data_utils.save_var(pseudo_obs,         file=os.path.join(self.snpe_folder, 'pseudo_obs.pkl'))
        data_utils.save_var(kernel_bandwidths,  file=os.path.join(self.snpe_folder, 'kernel_bandwidths.pkl'))
      
      # Make sure last posterior does exist.
      assert sample_distributions[-1] is not None
      
      print('\n----------------------------------------------------\n')
      
      i_round += 1

    print('---> Done!')
    
  ############################################################################
  def run_SNPE_round(
      self, n_train, proposal=None, verbose=True, initial_tds=None,
    ):
    ''' Run a single rounds of SNPE.
    
    Parameters:
    
    n_train : int
      Number of samples for this round.
      
    proposal : distribution or None
      Sample distribution or None for SNPE-C
      
    verbose : bool
      Print user information
      
    initial_tds : None or training data
      If round was started before and is now started again, use these samples.
      
    Returns:
    
    log, td, posterior : output of SNPE
    
    '''
    
    from delfi import __version__ as delfi_version
    
    # Update generator.
    self.inf_snpe.generator = self.create_delfi_generator()
    
    print(f'Inference round, using SNPE-{self.snpe_type}')
    
    # SNPE-B, modified.
    if self.snpe_type in ['b', 'B']:
      if delfi_version == '0.5.1':
        log_list, td_list, posterior_list = self.inf_snpe.run(
          n_train=n_train,
          n_rounds=1,
          epochs=self.nn_epochs,
          minibatch=self.nn_minibatch,
          round_cl=self.round_cl,
          stop_on_nan=True,
          proposal=proposal,
          initial_tds=initial_tds,
          verbose=verbose,
        )
      else:
        log_list, td_list, posterior_list = self.inf_snpe.run(
          n_train=n_train,
          n_rounds=1,
          epochs=self.nn_epochs,
          minibatch=self.nn_minibatch,
          round_cl=self.round_cl,
          stop_on_nan=True,
          proposal=proposal,
          verbose=False,
        )

    # SNPE-C.
    elif self.snpe_type in ['c', 'C']:
      assert delfi_version != '0.5.1', 'Use the latest version if you want to use this function.'
      
      # Run
      log_list, td_list, posterior_list = self.inf_snpe.run(
        n_train=n_train,
        n_rounds=1,
        epochs=self.nn_epochs,
        minibatch=self.nn_minibatch,
        round_cl=self.round_cl,
        stop_on_nan=True,
        proposal=self.snpec_proposal,
      )
    
    return log_list[0], td_list[0], posterior_list[0]
    
  ############################################################################
  def update_samples_file(self, samples_file):
    ''' Update the filename in which samples are stored.
    '''
    
    if samples_file is None:
      self.samples_file = os.path.join(self.samples_folder, 'delfi_samples.pkl')
    elif 'optim_data' in samples_file:
      self.samples_file = samples_file
    else:
      self.samples_file = os.path.join(self.samples_folder, samples_file)
    
    if self.samples_file[-4:] != '.pkl': self.samples_file += '.pkl'
    
  ############################################################################
  def load_tds_from_file(self, file, params):
    ''' Load training data to pass it to SNPE.
    Importance weights will be recomputed.
    '''
  
    print(file)
    loaded_samples = data_utils.load_var(file)
    
    if 'wall-time' in loaded_samples.keys():
      n_loaded_samples = loaded_samples['wall-time'].size
    elif 'loss' in loaded_samples.keys():
      n_loaded_samples = loaded_samples['loss']['total'].size
    else:
      raise NotImplementedError
      
    print('Loaded', n_loaded_samples, 'samples from', file)
    
    loss_names = loaded_samples['loss'].keys()
    tds_loss = []
    for i in range(n_loaded_samples):
      sample_loss_dict = {}
      for loss_name in loss_names:
        sample_loss_dict[loss_name] = loaded_samples['loss'][loss_name][i]
      tds_loss.append(self.to_network_input(sample_loss_dict))
    tds_loss = np.array(tds_loss)

    if tds_loss.ndim == 1: tds_loss = np.atleast_2d(tds_loss).T

    if isinstance(loaded_samples['params'], dict):
      n_params = len(loaded_samples['params'])
    else:
      n_params = loaded_samples['params'].shape[1]
    
    assert n_params == params.p_N
    
    tds_params = np.zeros((n_loaded_samples, n_params))
    
    if isinstance(loaded_samples['params'], dict):
      for idx, param in enumerate(params.p_names):
        tds_params[:,idx] = params.sim_param2opt_param(loaded_samples['params'][param], param)
    else:
      tds_params = loaded_samples['params']
      
    loaded_tds = (tds_params, tds_loss, None)
    return loaded_tds, n_loaded_samples
    
  ############################################################################
  def load_SNPE_rounds(self):
    ''' Load stored rounds of SNPE.
    '''

    inf_snpes            = data_utils.load_var(os.path.join(self.snpe_folder, 'inf_snpes.pkl'))
    logs                 = data_utils.load_var(os.path.join(self.snpe_folder, 'logs.pkl'))
    tds                  = data_utils.load_var(os.path.join(self.snpe_folder, 'tds.pkl'))
    sample_distributions = data_utils.load_var(os.path.join(self.snpe_folder, 'sample_distributions.pkl'))
    n_samples            = data_utils.load_var(os.path.join(self.snpe_folder, 'n_samples.pkl'))
      
    if self.snpe_type in ['b', 'B']:
      kernel_bandwidths = data_utils.load_var(os.path.join(self.snpe_folder, 'kernel_bandwidths.pkl'))
      pseudo_obs        = data_utils.load_var(os.path.join(self.snpe_folder, 'pseudo_obs.pkl'))
    else:
      kernel_bandwidths = None
      pseudo_obs = None
    
    print('Loaded {:d} rounds'.format(len(n_samples)))
    
    return inf_snpes, logs, tds, sample_distributions, n_samples, kernel_bandwidths, pseudo_obs
    
  ############################################################################
  def create_SNPE_rounds_backup(self):
    ''' Moves all data of previous rounds to a backup folder.
    '''
  
    # Move old sample files to not overwrite them.
    samples_files = os.listdir(self.samples_folder)
    samples_files = [os.path.join(self.samples_folder, samples_file) for samples_file in samples_files]
    for samples_file in samples_files:
      backup_samples_file = samples_file.replace('samples/', 'backups/')
      # Delete old backup file if necessary. Store maximally one backup file.
      if os.path.isfile(backup_samples_file): os.remove(backup_samples_file)
      os.rename(samples_file, backup_samples_file)
      
    # Move old snpe files.
    snpe_files = os.listdir(self.snpe_folder)
    snpe_files = [self.snpe_folder + snpe_file for snpe_file in snpe_files]
    for snpe_file in snpe_files:
      backup_snpe_file = snpe_file.replace('snpe/', 'backups/')
      # Delete old backup file if necessary. Store maximally one backup file.
      if os.path.isfile(backup_snpe_file): os.remove(backup_snpe_file)
      os.rename(snpe_file, backup_snpe_file)
  
  ############################################################################
  def gen_samples(
      self, params_list, n_reps=1, verbose=True, save_data=True, **kwargs
    ):
    ''' Generate samples for given parameters.
    '''
    assert n_reps == 1, 'n_reps > 1 not implemented yet'
    
    if verbose: print('\t\tDrawing ' + str(len(params_list)) + ' samples ... ', end='')
    self.load_random_state()
  
    model_output_list = self.optim.run_parallel(params_list, save_data=False, **kwargs)
    if verbose: self.print_total_loss_info(model_output_list)
    
    if save_data:
      self.optim.save_samples(model_output_list, file=self.samples_file, overwrite=False)
      self.save_random_state()
    else:
      if verbose: print(f'Not saving {len(params_list)} samples.')
    
    network_samples = [[self.to_network_input(model_output['loss'])] for model_output in model_output_list]
    return network_samples
    
  ############################################################################
  def to_network_input(self, model_output_loss):
    ''' Transform data to training data for SNPE.
    '''
    if self.scalar_loss:
      model_output_loss_for_nn = model_output_loss['total']
    else:
      model_output_loss_for_nn = np.array([v for k, v in model_output_loss.items() if k != 'total'])
      
    return model_output_loss_for_nn
    
  ############################################################################
  @staticmethod
  def print_total_loss_info(model_output_list):
    ''' Print min, mean and max totoal loss.
    '''
    loss_arr = np.array([model_output['loss']['total'] for model_output in model_output_list])
    min_loss  = np.nanmin(loss_arr)
    mean_loss = np.nanmean(loss_arr)
    max_loss  = np.nanmax(loss_arr)
    print('Total-Loss: Min={:.3f}; Mean={:.3f}; Max={:.3f}; NaN={:.2f}%'.format(\
      min_loss, mean_loss, max_loss, 100.*np.sum(np.isnan(loss_arr))/loss_arr.size
    ))
    
  ############################################################################
  def load_samples(
      self, files=None, verbose=False, concat_traces=True, list_traces=False,
      return_n_samples=False, return_sort_idx=False,
    ):
    ''' Loads samples and stacks them. Can also list them is sizes differ.
    TODO: Make cleaner.
    '''

    assert not (concat_traces and list_traces), 'choose one'
  
    # If no files are given, use default.
    if files is None: files = os.listdir(self.samples_folder)
    
    # Append folder.
    files = sorted([os.path.join(self.samples_folder, file) if 'optim_data' not in file else file for file in files])
    
    # Load all files.
    loaded_samples_list = [self.optim.load_samples(file=file, verbose=verbose) for file in files]
    
    if return_n_samples:
      n_samples = np.cumsum([loaded_samples['loss']['total'].shape[0] for loaded_samples in loaded_samples_list if loaded_samples is not None])    
      if verbose: print('Loaded ' + str(n_samples[-1]) + ' samples!')

    # Use first file to initialize.
    concatenated_samples = loaded_samples_list[0]
    params               = list(concatenated_samples.keys())
    
    # Handle special parameters.
    for param in params:
      if param not in ['loss', 'params', 'wall-time']:
        if list_traces:
          if not isinstance(concatenated_samples[param], list):
            concatenated_samples[param] = [concatenated_samples[param]]
        elif not concat_traces:
          # Remove.
          del concatenated_samples[param]
    
    # Append all other files.
    for loaded_samples in loaded_samples_list[1:]:
      if loaded_samples is not None:
        
        # Append loss.
        for loss_name in concatenated_samples['loss'].keys():
          if loss_name in loaded_samples['loss'].keys():
            concatenated_samples['loss'][loss_name] = np.concatenate((concatenated_samples['loss'][loss_name], loaded_samples['loss'][loss_name]))
        
        # Append params.
        if isinstance(concatenated_samples['params'], dict):
          for opt_param in concatenated_samples['params'].keys():
            if opt_param in loaded_samples['params'].keys():
              concatenated_samples['params'][opt_param] = np.concatenate((concatenated_samples['params'][opt_param], loaded_samples['params'][opt_param]))
        elif isinstance(concatenated_samples['params'], np.ndarray):
          concatenated_samples['params'] = np.concatenate((concatenated_samples['params'], loaded_samples['params']))

        # Append wall-time.
        if 'wall-time' in concatenated_samples:
          concatenated_samples['wall-time'] = np.concatenate((concatenated_samples['wall-time'], loaded_samples['wall-time']))

        # Concatenate rest.
        for param in concatenated_samples:
          if param not in ['loss', 'params', 'wall-time']:
            if concat_traces:
              if isinstance(concatenated_samples[param], dict):
                for param2 in concatenated_samples[param].keys():
                  if param2 in loaded_samples[param].keys():
                    concatenated_samples[param][param2] = np.concatenate((concatenated_samples[param][param2], loaded_samples[param][param2]))
              elif isinstance(concatenated_samples[param], np.ndarray):
                concatenated_samples[param] = np.concatenate((concatenated_samples[param], loaded_samples[param]))
            elif list_traces:
              if isinstance(loaded_samples[param], list):
                concatenated_samples[param] += loaded_samples[param]
              else:
                concatenated_samples[param] += [loaded_samples[param]]
    
    # Return.
    return_values = [concatenated_samples]
    if return_n_samples: return_values += [n_samples]
    if return_sort_idx:  return_values += [np.argsort(concatenated_samples['loss']['total'])]
    return return_values
    
############################################################################
############################################################################
class EmptyDELFI_Optimizer(DELFI_Optimizer):

  def __init__(self, *args, **kwargs):
    self.optim = optim_funcs.EmptyOptimizer()
    pass
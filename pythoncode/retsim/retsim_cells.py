import pandas as pd
import numpy as np
import warnings
import os
import copy
from multiprocessing import Pool
from matplotlib import pyplot as plt

import python2retsim
import plot_cell_rec_data
import interpolation_utils
import data_utils

XCONE  = 0
DBP1   = 6
HBP1   = 10
GCA    = 23
GCAOFF = 26
AII    = 13

############################################################################
############################################################################
class Cell():
  
  ############################################################################
  def __init__(
      self, timeout=3600,
      retsim_path='/gpfs01/berens/user/joesterle/berens/neuronc-master/nc/models/retsim/'
    ):
    
    # Input file names.
    self.expt_file_list = None
   
    self.bp_morphfile_ON  = None
    self.bp_morphfile_OFF = None
    self.bp_densfile_ON = None
    self.bp_densfile_OFF = None
    self.cone_densfile = None
    
    self.nval_file = None
    self.chanparams_file = None
    self.bp_type = None
    
    # Output file names.
    self.compfile = 'temp_comp_file'
    self.comsol_compfile = None
    
    # Morphology parameters and data.
    self.mxrot = -90
    self.myrot = 0
    self.mzrot = 180
    self.comsol_x0 = 0
    self.comsol_y0 = 0
    self.comsol_z0 = 0
    self.cone_soma_z = 30
    self.morph_data = None
    self.comp_data = None
    self.auto_correct_comp_data = True
    self.n_bc_comps = None
    self.n_gc_comps = None
    self.n_cones = None
    
    # Other simulation parameters.
    self.predur = 2
    self.rec_dt = 1e-3
    self.sim_dt = 1e-5
    self.predur_dt = 1e-4
    self.stim_dt = 1e-4
    self.syn_dt = 1e-4
    self.timeout = timeout
     
    # Stimulus data.
    self.stim = None
    self.retsim_stim_file_base = "temp_stim"
    self.retsim_stim_file = None
    
    # Stimulus data in t_rng.
    self.__t_rng = None
    self.__t_rng_from_user = False
    self.stim_in_t_rng = None
    
    # Recorded data.
    self.rec_data = {}

    # Retsim output image.
    self.retsim_im = None
    
    # Simulation parameters.
    self.params_default = {}
    self.params_unit = {}
    
    self.is_OFF_bp = False
    self.make_dbp1 = True
    self.make_hbp1 = False
    
    self.set_n_cones = 0
    self.make_cones = False
    self.merge_cones = True
    
    self.set_retsim_path(retsim_path)
    
    os.path.join(self.retsim_path, 'comp_files')
    os.path.join(self.retsim_path, 'stim_files')

  ############################################################################
  def set_retsim_path(self, retsim_path):
    ''' Set retsim path, i.e. tell python where retsim is installed.
    '''
    self.retsim_path = retsim_path
    python2retsim.set_dir(self.retsim_path)
  
  ############################################################################  
  def set_stim(self, stim, stim_type=None, verbose=False, n_reps=1, inter_rep_dt=0.0):
    ''' Set stimulus or stimuli.
    
    Paramters:
    
    stim : DataFrame or list of DataFrames
      Stimulus for cell.
      
    stim_type : str
      Stimulus types, e.g. 'Light' or 'Vext'
      Note that different stimulus types need different types of DataFrames.
      
    verbose : bool
      Print information
      
    n_reps : int
      Repeat stimulus, default to 1, i.e. no repetition.
      
    inter_rep_dt : float
      Time delay between two repetitions.
    '''
    if stim is None:
      self.stim = None
      return
  
    # Make stimulus to list.
    if not(isinstance(stim, list)): stim = [stim]
    
    # Get stimulus type
    stim_type = stim_type or self.stim_type
    assert stim_type in ['Light', 'Vclamp', 'Iclamp', 'Vext'], 'Not implemented'
    self.stim_type = stim_type
    
    for stim_i in stim:
      self.check_stimulus_format(stim_i)
    
    if n_reps == 1:
      self.stim = [stim_i.astype(float) for stim_i in stim]
    
    else:
      self.stim = []
      for stim_i in stim:
        self.stim.append(self.repeat_stimulus(stim_i, n_reps=n_reps, inter_rep_dt=inter_rep_dt))
    
    # Update t_rng, stimulus duration and stimulus in t_rng.
    self.update_t_rng() 
  
  ############################################################################
  def check_stimulus_format(self, stim_df):
    ''' Test if stimulus has correct format.
    '''
    assert isinstance(stim_df, pd.DataFrame), 'Stimuli must be DataFrames'
    
    if self.stim_type in ['Light', 'Vclamp', 'Iclamp']:
      assert set(stim_df.columns) == {'Time', 'Stim'}
    
    elif self.stim_type == 'Vext':
      assert stim_df.columns[0] == 'Time', 'Fist columns must be Time'
      if self.n_bc_comps is not None: 
        if (stim_df.shape[1]-1 != self.n_bc_comps) and (stim_df.shape[1]-1 != self.n_bc_comps + self.n_cone_comps):
          warnings.warn('{:d} External compartments vs. '
                        '{:d} BC compartments and '
                        '{:d} cone compartments'.format(
                        stim_df.shape[1]-1, self.n_bc_comps, self.n_cone_comps))
  
  ############################################################################
  @staticmethod
  def repeat_stimulus(self, stim_df, n_reps, inter_rep_dt):
    ''' Repeat a stimulus.
    
    Parameters:
    
    stim_df : DataFrame
      Stimulus DataFrame
    
    n_reps : int
      Number of repetitions.
      
    inter_rep_dt : float
      Time delay between two repetitions.
    
    Returns:
    
    stim_df : DataFrame
      Repeated stimulus DataFrame
    
    '''
    n_rows = stim_df.shape[0]
      
    # Repeat stimulus.
    stim_df = stim_df.append([stim_df]*(n_reps-1), ignore_index=True)

    # Get time difference between repetitions.
    rep_dt = stim_df['Time'][n_rows-1] + stim_df['Time'][1] + inter_rep_dt
    
    # Adjust time.
    for i in range(n_reps-1):
      stim_df['Time'].iloc[(i+1)*n_rows:] += rep_dt

    return stim_df
  
  ############################################################################
  def update_t_rng(self, t_rng=None):
    ''' Set a new range, i.e. a starting and stopping point for the simulation.
    If necessary, stimuli are updated.
    '''
    
    # Update t_rng.
    if t_rng is not None:
      self.__t_rng = t_rng
      self.__t_rng_from_user = True
    elif (self.stim is not None) and not(self.__t_rng_from_user):
      self.__t_rng = (0, np.max(self.stim[0]['Time']))
    
    if self.__t_rng is not None:
      # Update stimulus in t_rng.
      if self.stim is not None:
        self.stim_in_t_rng = [self.__data_frame2t_rng(stim_i, self.__t_rng)[0] for stim_i in self.stim]
      
      # Update stimulus duration.
      self.stimdur = self.__t_rng[1] - self.__t_rng[0]
      
      # Update retsim stimulus file.
      self.retsim_stim_file = self.retsim_stim_file_base + '_{:.2f}-{:.2f}'.format(self.__t_rng[0], self.__t_rng[1])
    
    return self.__t_rng
  
  ############################################################################
  def update_cpl(self, cpl=None, verbose=True, **kwargs):
    ''' Update cpl parameter of cell.
    cpl sets the number of compartments.
    Updates retsim object.
    
    Parameters:
    
    cpl : float
      New cpl parameter.
      
    kwargs
      Additional cpl parameters.
      
    verbose : bool
      If True, will print user information.
    
    '''
    
    must_recompile = False
    
    # Use cpl for all cpl values.
    if cpl is not None:
      assert isinstance(cpl, (float, int))
      found_cpl = False
      for k,v in self.params_default.items():
        if k[0:3] == 'cpl':
          self.params_default[k] = cpl
          found_cpl = True
      if not found_cpl:
        print('No cpl param found, will use default cpl_a and cpl_d.')
        must_recompile = True
        self.params_default['cpl_a'] = cpl
        self.params_default['cpl_d'] = cpl
    
    # Use specific cpl values. Not alternative, might use also cpl.
    if len(kwargs) > 0:
      for k,v in kwargs.items():
        if k[0:3] == 'cpl':
          if k not in self.params_default.keys():
            print('Added new CPL parameter ' + k + ' you need to recompile expt file')
            must_recompile = True
          self.params_default[k] = kwargs[k]
    self.init_retsim(verbose=False, update=True)
    
    if verbose: print('# New number of BC compartments: ' + str(self.n_bc_comps))
  
  ############################################################################
  @staticmethod
  def __data_frame2t_rng(data_frame, t_rng):
    ''' Fit data frame to t_rng.
    '''
    # Get time course of stimulus.
    time = data_frame['Time']
    
    new_time = np.append(t_rng[0], np.append(time[np.logical_and(time > t_rng[0], time < t_rng[1])], t_rng[1]))   
    new_data_frame = pd.DataFrame({'Time': new_time - t_rng[0]})
    
    for column in data_frame.columns:
      if column != 'Time':
        new_data_frame[column] = interpolation_utils.in_ex_polate(x_old=time, y_old=data_frame[column], x_new=new_time)

    return new_data_frame, t_rng   

  ############################################################################
  def create_retsim_stim_file_name(self, stim_idx=None):
    ''' Create a file name for a retsim stimulus.
    '''
    
    if stim_idx is None: suffix = '.csv'
    else:                suffix = '_' + str(stim_idx) + '.csv'
    return os.path.join(self.retsim_path, 'stim_files', self.retsim_stim_file + suffix)
    
  ############################################################################
  def create_retsim_stim_file(self, t0=0, stim_idx=None, verbose=False):
    ''' Create a retsim stimulus file.
    '''
    
    if stim_idx is None:
      assert len(self.stim) == 1, 'For multiple stimuli, please state stimulus index'
      stim_i = self.stim_in_t_rng[0].copy()
    else:
      stim_i = self.stim_in_t_rng[stim_idx].copy()
    
    output_file = self.create_retsim_stim_file_name(stim_idx)
    if verbose: print(output_file)
    
    if self.stim_type == 'Vext':
      self.__create_retsim_stim_file_Vext(stim=stim_i, output_file=output_file, t0=t0, verbose=verbose)
    else:
      self.__create_retsim_stim_file_other(stim=stim_i, output_file=output_file, t0=t0)
    
  
  ############################################################################
  @staticmethod
  def __create_retsim_stim_file_other(stim, output_file, t0):
    ''' Create a retsim stimulus file, Light etc.
    '''
    stim['Time'] += t0
    stim.to_csv(output_file, index=False, header=False, sep=',')
  
  ############################################################################  
  def __create_retsim_stim_file_Vext(
      self, stim, output_file, t0, verbose=False, inc_cones=True, inc_bc=True, inc_gc=False,
    ):
    ''' Create a retsim stimulus file, extracellular voltage.
    '''
    
    stim = stim.copy()
    
    assert stim.columns[0] == 'Time', 'first column must be Time followed by Vext for every comp'
   
    # Get compartment indices.
    valid_idxs = np.zeros(self.comp_data.shape[0], dtype=bool)
    
    cone_idxs = (self.comp_data['cell_type'] == XCONE)
    bp_idxs   = (self.comp_data['cell_type'] == DBP1) | (self.comp_data['cell_type'] == HBP1)
    gc_idxs   = (self.comp_data['cell_type'] == GCA) | (self.comp_data['cell_type'] == GCAOFF)
    
    if inc_cones: valid_idxs = valid_idxs | cone_idxs
    if inc_bc:    valid_idxs = valid_idxs | bp_idxs
    if inc_gc:    valid_idxs = valid_idxs | gc_idxs
   
    # Get number of compartments according to stimulus and stimulus length.
    n_comps = np.sum(valid_idxs)
    assert n_comps == (stim.shape[1] - 1)
    
    if verbose: print("Vext for " + str(n_comps) + " compartments" , end='')
    n_time = stim.shape[0]
    
    # Generate time vector.
    ncstim_t = t0 + np.repeat(stim['Time'].values, 2*n_comps)
    
    # Generate mask vector. Alternating zeros and ones.
    ncstim_mask = np.tile(np.append(np.zeros(n_comps, dtype=int), np.ones(n_comps, dtype=int)), n_time)
    
    # Generate voltage vector.
    ncstim_v = np.tile(stim.iloc[:,1:].values, 2).flatten()
    ncstim_v *= np.tile(np.append(np.ones(n_comps, dtype=int), -1*np.ones(n_comps, dtype=int)), n_time)
    
    # Append end.
    ncstim_t    = np.append(ncstim_t, (t0 + stim['Time'].values[-1] + self.rec_dt) * np.ones((2*n_comps, 1)))
    ncstim_v    = np.append(ncstim_v, np.zeros((n_comps, 1), dtype=int))
    ncstim_mask = np.append(ncstim_mask, np.zeros((n_comps, 1), dtype=int))
    
    # Append start
    if t0 > 0:
      ncstim_t    = np.append(np.zeros((n_comps, 1)), ncstim_t)
      ncstim_v    = np.append(np.zeros((2*n_comps, 1)), ncstim_v)
      ncstim_mask = np.append(np.append(np.zeros(n_comps, dtype=int), np.ones(n_comps, dtype=int)), ncstim_mask)
    else:
      ncstim_t = ncstim_t[n_comps:]
    
    # Generate rest.
    n_reps = int(ncstim_t.shape[0] / n_comps)
    ncstim_nd1 = np.tile(self.comp_data['cell_type'][valid_idxs].values, n_reps)
    ncstim_nd2 = np.tile(self.comp_data['cell_number'][valid_idxs].values, n_reps)      
    ncstim_nd3 = np.tile(np.arange(5000, 5000+n_comps), n_reps)
    
    ncstim_int = -1 * np.ones(ncstim_t.shape, dtype=int)
    ncstim_act = np.chararray(ncstim_t.shape, itemsize=7)
    ncstim_act[:] = 'evclamp'
    ncstim_sequence = np.zeros(ncstim_t.shape, dtype=int)

    # Create data frame.
    stim_neuronc = pd.DataFrame(
      {
        'time': ncstim_t,
        'ndim1': ncstim_nd1,
        'ndim2': ncstim_nd2,
        'ndim3': ncstim_nd3,
        'intens': ncstim_int,
        'wavel': ncstim_v,
        'mask': ncstim_mask,
        'action': ncstim_act.astype('str'),
        'sequence': ncstim_sequence
       },
       columns = ['time', 'ndim1', 'ndim2', 'ndim3', 'intens', 'wavel', 'mask', 'action', 'sequence']
    )
    
    # Write to file.
    if verbose: print("Save stimulus to\n" + output_file)
    stim_neuronc.to_csv(
      output_file, sep = '\t', float_format='%.9f', index = False, \
      header = ["#time", "node", "", "", "intens", "wavel", "mask", "action", "sequence"]
    )
  
  ############################################################################
  def create_retsim_expt_file(
      self, base_file=None, output_file=None, to_file=True,
      set_n_cones=None, set_n_cones_on=None, set_n_cones_off=None, verbose=True,
      deactivate_outer_segment=False, on2cone_nodes=None, off2cone_nodes=None,
      off_cell=None, on_cell=None, dxdy_on=[0,0], dxdy_off=[0,0],
    ):
    
    ''' Create a retsim experiment file.
    '''
    
    # Get base file.
    if base_file is not None:
      if isinstance(base_file, list): base_file_list = base_file
      else:                           base_file_list = [base_file]
    else:
      base_file_list = self.expt_base_file_list
    
    if verbose: print('Base files: ' + str(base_file_list))
    
    # Get output file. (= The file were the experiment will be stored)
    if output_file is not None:
      if isinstance(output_file, list): outputfile_list = output_file
      else:                             outputfile_list = [output_file]
    else:
      outputfile_list = self.expt_file_list

    # Assert #input is equal to #output.
    assert len(base_file_list) == len(outputfile_list), 'must be equal'
    
    with open(os.path.join(self.retsim_path, 'makefile'), 'r') as makefile:
      makefiledata = makefile.read()
    
    for file in outputfile_list:
      if f"expt_{file}.$(DLSUFFIX)" not in makefiledata:
        print(f'Warning: Experiment {file} is not in retsim makefile. '
              'Therefore it will not be compiled per default. '
              'This might causes errors.')
    
    # Create all output files.
    for idx in range(len(base_file_list)):
      base_file = base_file_list[idx]
      output_file = outputfile_list[idx]
      if set_n_cones_on is None:
        if set_n_cones is None:
          set_n_cones_on = self.set_n_cones
        else:
          set_n_cones_on = set_n_cones
      if set_n_cones_off is None:
        if set_n_cones is None:
          set_n_cones_off = self.set_n_cones
        else:
          set_n_cones_off = set_n_cones
      
      with open(base_file, 'r') as f:
        base_lines = f.readlines()
      
      # Find lines in base experiment file to paste parameters.
      found_declare = False
      found_init = False
      found_print = False
      for idx, line in enumerate(base_lines):
        if line == 'declare_params_here\n':
          found_declare = True
          idx_declare = idx
        if line == '  init_params_here\n':
          found_init = True
          idx_init = idx
        if line == '  // Add printer command here.\n':
          found_print = True
          idx_print = idx
            
      assert found_declare and found_init
      
      # Create new lines.
      skip_params_declare = [
        "n_c_bp", "rec_predur", "post_load_predur_timinc",
        "load_eq", "save_eq", "run_predur_only", "eqfile",
        "merge_cones", "set_dscaeg", "set_c_vst", "set_tempcel",
        "post_load_predur", "pre_load_predur",
        "set_cone_timec", "set_cone_loopg",
      ] # Skip declaring this params.
      
      params_declare_list = []
      for param, value in self.params_default.items():
        if param not in skip_params_declare:
          if isinstance(value,(int, float)):
            params_declare_list.append('double ' + param + ';\n')
          elif isinstance(value, dict) and 'cpl' in param:
            params_declare_list.append('double ' + param + ';\n')
      
      params_init_list = []
      for param, value in self.params_default.items():
        if param not in skip_params_declare and isinstance(value,(int, float)):
          params_init_list.append(('  setptr(\"'+ param + '\",').ljust(25) + ('&' + param + ');').ljust(20)\
                                  + param.ljust(15) + ' = {:.4g};\n'.format(value * self.get_unit(param)))

      print_lines = []
      if found_print:
        for param in self.params_default.keys():
          print_lines.append('  PRINTER(' + param + '); fprintf(stderr, "\t");\n')
      else:
        idx_print = -1
    
    
      # Concatenate lines.
      outputlines = base_lines[0:idx_declare] +\
                    params_declare_list +\
                    base_lines[idx_declare+1:idx_init] +\
                    params_init_list +\
                    base_lines[idx_init+1:idx_print] +\
                    print_lines +\
                    base_lines[idx_print:]
      
      # Get cone connections and cone positions.
      if on2cone_nodes: # Neither None nor empty list.
        assert isinstance(on2cone_nodes, (list, np.ndarray))
      elif (on_cell is not None) and self.make_cones:
        on2cone_nodes = on_cell.get_bc2cone_nodes(verbose=False, set_n_cones=set_n_cones_off)
      elif (self.bp_morphfile_ON is not None) and self.make_cones:
        on2cone_nodes = self.get_bc2cone_nodes(verbose=False, set_n_cones=set_n_cones_on, morphfile=self.bp_morphfile_ON, bp_type='BC5')
      else:
        on2cone_nodes = []
            
      if off2cone_nodes: # Neither None nor empty list.
        assert isinstance(on2cone_nodes, (list, np.ndarray))
      elif (off_cell is not None) and self.make_cones:
        off2cone_nodes = off_cell.get_bc2cone_nodes(verbose=False, set_n_cones=set_n_cones_off)
      elif (self.bp_morphfile_OFF is not None) and self.make_cones:
        off2cone_nodes = self.get_bc2cone_nodes(verbose=False, set_n_cones=set_n_cones_off, morphfile=self.bp_morphfile_OFF, bp_type='BC3')
      else:
        off2cone_nodes = []
      
      # Get node positions.
      if np.asarray(on2cone_nodes).size > 0:
        if on_cell is not None:
          on2cone_nodes_x,  on2cone_nodes_y, _ = on_cell.get_bc_node_retsim_xyz(on2cone_nodes)
        else:
          on2cone_nodes_x,  on2cone_nodes_y, _ = self.get_bc_node_retsim_xyz(on2cone_nodes, morphfile=self.bp_morphfile_ON)
          
        on2cone_nodes_x += dxdy_on[0]
        on2cone_nodes_y += dxdy_on[1]
      else:
        on2cone_nodes_x, on2cone_nodes_y = [], []
        
      if np.asarray(off2cone_nodes).size > 0:
        if off_cell is not None:
          off2cone_nodes_x, off2cone_nodes_y, _ = off_cell.get_bc_node_retsim_xyz(off2cone_nodes)
        else:
          off2cone_nodes_x, off2cone_nodes_y, _ = self.get_bc_node_retsim_xyz(off2cone_nodes, morphfile=self.bp_morphfile_OFF)
        
        off2cone_nodes_x += dxdy_off[0]
        off2cone_nodes_y += dxdy_off[1]
      else:
        off2cone_nodes_x, off2cone_nodes_y = [], []
      
      # Lines to replace in experiment file.
      lines2replace = {
        "int contact_nodes_dbp1 [] = ": "int contact_nodes_dbp1 [] = {" + str(list(on2cone_nodes))[1:-1] + "};\n",
        "int contact_nodes_hbp1 [] = ": "int contact_nodes_hbp1 [] = {" + str(list(off2cone_nodes))[1:-1] + "};\n",
        
        "double cone_x_dbp1 []": "double cone_x_dbp1 [] = {" + str(list(on2cone_nodes_x))[1:-1] + "};\n",
        "double cone_y_dbp1 []": "double cone_y_dbp1 [] = {" + str(list(on2cone_nodes_y))[1:-1] + "};\n",
        
        "double cone_x_hbp1 []": "double cone_x_hbp1 [] = {" + str(list(off2cone_nodes_x))[1:-1] + "};\n",
        "double cone_y_hbp1 []": "double cone_y_hbp1 [] = {" + str(list(off2cone_nodes_y))[1:-1] + "};\n",
      }
        
      
      # Deactivate output segment by setting conductance to zero.
      if deactivate_outer_segment:
        lines2replace['cone_maxcond'] = "  cone_maxcond = 0.0;\n"
      
      # Replace.
      for idx, line in enumerate(outputlines):
        for key, value in lines2replace.items():
          if key in line:
            outputlines[idx] = value
            break
      
      if to_file:
        # Write to file.
        with open(self.retsim_path + 'expt_' + output_file + '.cc', 'w') as f:
          for line in outputlines:
              f.write(line)
      else:
        return outputlines
        
    if verbose: print('Output files: ' + str(outputfile_list))

  ############################################################################  
  def compute_area(self, regions=None):
    ''' Compute the membrane surface area of the cell.
    
    Parameters:
    
    regions : list of str or None
      Regions to consider.
      If None, will consider all regions. Default is None.
      
    Returns:
      
    area : float
      Area of the membrane in cm^2
      
    '''
  
    if regions is None:
      regions = np.unique(self.morph_data['region'])

    A = 0
    for node_idx in range(self.morph_data.shape[0]):
        if self.morph_data['region'][node_idx] in regions:
            x1 = self.morph_data['x'][node_idx]
            y1 = self.morph_data['y'][node_idx]
            z1 = self.morph_data['z'][node_idx]
            for p_idx in self.morph_data['p_idx'][node_idx]:
                x2 = self.morph_data['x'][p_idx]
                y2 = self.morph_data['y'][p_idx]
                z2 = self.morph_data['z'][p_idx]
    
            A += np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) * np.pi * self.morph_data['dia'][node_idx]
    if "SOMA" in regions:
      A += (self.morph_data['dia'][self.morph_data['region'] == "SOMA"][0]/2)**2 * np.pi * 4
    print("A = {:.2f} um²".format(A))
    print("A = {:.2g} cm²".format(A*1e-8))
    return A*1e-8
  
  ############################################################################
  def read_retsim_morphfile(self, morphfile=None, verbose=False, update=True):
    ''' Read the morphfile created by retsim.
    
    Parameters:
    
    morphfile : str
      Name of the morphfile
      
    verbose : bool
      If True, will print user information.
      
    update : bool
      If True, will update cell object with read data.
    
    '''

    assert self.retsim_path != None
    
    if morphfile is None:
      assert (self.bp_morphfile_ON is None) != (self.bp_morphfile_OFF is None) # XOR
      morphfile = self.bp_morphfile_ON or self.bp_morphfile_OFF
  
    # Read file.
    morph_data = pd.read_csv(
      self.retsim_path + 'runconf/' + morphfile, delim_whitespace=True, header=None,
      comment='#', names=['node','parent','dia','x','y','z','region','dend']
    )
    
    # Load morphology file.
    n_nodes = morph_data.shape[0]    
    
    # Print.
    if verbose:
      print('Raw morph_data:')
      for dim in ['x', 'y', 'z']: print(dim + ' in ({:.2f}, {:.2f})'.format(np.min(morph_data[dim]), np.max(morph_data[dim])), end='\t')
      print()

    # For every node: Find parent nodes. Every node must have a parent node which can be itself.
    morph_data['p_idx'] = np.empty((n_nodes, 0)).tolist()
    for i in range(n_nodes):
      j = np.argwhere(morph_data['parent'].values[i] == morph_data['node'].values)
      morph_data.loc[i, 'p_idx'] = list(j)
    
    # For every node: Find child nodes. If there is no child, point to node itself.
    morph_data['c_idx'] = np.empty((n_nodes, 0)).tolist()
    for i in range(n_nodes):
      j = np.argwhere(morph_data['node'].values[i] == morph_data['parent'].values)
      if j.size > 0: morph_data.loc[i, 'c_idx'] = list(j)
      else:          morph_data.loc[i, 'c_idx'] = [[i]]
    
    if update: self.morph_data = morph_data
    
    return morph_data
   
  ############################################################################
  def add_comsol_xyz_to_morph_data(self, verbose=True):  
  
    x = self.morph_data['x'].values
    y = self.morph_data['y'].values
    z = self.morph_data['z'].values
  
    if verbose:
      print('\nRaw morphology data:')
      print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
      print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
      print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)))
    
    # Rotate.
    x, y, z = self.rotate_xyz(x, y, z)
    
    # Add offset.
    x, y, z = self.offset_xyz(x, y, z)
    
    if verbose:
      print('\nCOMSOL morphology data.')
      print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
      print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
      print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)))

    # Save.
    self.morph_data['comsol_x'] = x
    self.morph_data['comsol_y'] = y  
    self.morph_data['comsol_z'] = z
  
  ############################################################################
  def read_retsim_compfile(self, compfile=None, verbose=True, update=True):
    
    # Load comps file.
    compfile = compfile or self.compfile
    
    comp_data_raw =  pd.read_csv(
      self.retsim_path + 'comp_files/' + compfile, delim_whitespace=True, comment='#', index_col=False,
      names=['is_comp','comp_idx', 'is_node', 'node_idx', 'nodes_to_comp', 'nodes',
      'is_xloc', 'x', 'is_yloc','y', 'is_zloc', 'z', 'is_dia', 'dia', 'is_Cm', 'Cm', 'is_Rm', 'Rm']
    )
   
    # Take only comp data values.
    comp_data = comp_data_raw[comp_data_raw['is_comp'] == 'comp']
   
    # Make lists out of strings.
    comp_nodes_list = []
    for row in range(comp_data.shape[0]):
      comp_nodes_raw = comp_data['nodes'].iloc[row][1:-1].split(';')[:-1]
      comp_nodes = []
      for comp_node in comp_nodes_raw:
        comp_nodes.append([int(i) for i in comp_node[1:-1].split(',')])
      if len(comp_nodes) == 0:
        comp_nodes = comp_nodes_list[-1]
      comp_nodes_list.append(comp_nodes)
      
    comp_data['nodes'] = [np.asarray(comp_nodes)[:,2] for comp_nodes in comp_nodes_list]
      
    # Assert only nodes from one cell per comp.
    cell_type = np.empty(comp_data.shape[0], dtype=int)
    cell_number = np.empty(comp_data.shape[0], dtype=int)
    
    no_nodes_idx_list = []
    
    for idx, comp_nodes in enumerate(comp_nodes_list):
      if len(comp_nodes) > 0:
        assert np.unique(np.asarray(comp_nodes)[:,0]).size == 1
        assert np.unique(np.asarray(comp_nodes)[:,1]).size == 1
        
        cell_type[idx] = comp_nodes[0][0]
        cell_number[idx] = comp_nodes[0][1]
      else:
        no_nodes_idx_list.append(idx)
    
    # Use previous.
    for idx in no_nodes_idx_list:
      cell_type[idx] = cell_type[idx-1]
      cell_number[idx] = cell_type[idx-1]
    
    comp_data['cell_type'] = cell_type
    comp_data['cell_number'] = cell_number

    comp_data = comp_data.drop(
      ['is_comp', 'is_node', 'nodes_to_comp', 'is_xloc','is_yloc', 'is_zloc', 'is_dia', 'is_Cm', 'is_Rm'],
      axis=1
    )
    comp_data = comp_data.set_index('comp_idx')
    
    # Count.
    n_cones        = np.unique(comp_data['cell_type'][comp_data['cell_type'] == XCONE]).size
    n_cone_comps   = np.sum(comp_data['cell_type'] == XCONE)
    n_dbp1_comps   = np.sum(comp_data['cell_type'] == DBP1)
    n_hbp1_comps   = np.sum(comp_data['cell_type'] == HBP1)
    n_gca_comps    = np.sum(comp_data['cell_type'] == GCA)
    n_gcaoff_comps = np.sum(comp_data['cell_type'] == GCAOFF)
    n_aii_comps    = np.sum(comp_data['cell_type'] == AII)
    
    # Print.
    if n_cones        and verbose > 0: print("n_cones        = " + str(n_cones))
    if n_cone_comps   and verbose > 0: print("n_cone_comps   = " + str(n_cone_comps))
    if n_dbp1_comps   and verbose > 0: print("n_dbp1_comps   = " + str(n_dbp1_comps))
    if n_hbp1_comps   and verbose > 0: print("n_hbp1_comps   = " + str(n_hbp1_comps))
    if n_gca_comps    and verbose > 0: print("n_gca_comps    = " + str(n_gca_comps))
    if n_gcaoff_comps and verbose > 0: print("n_gcaoff_comps = " + str(n_gcaoff_comps))
    if n_aii_comps    and verbose > 0: print("n_aii_comps    = " + str(n_aii_comps))
    
    # Update cell data?
    if update:
      self.n_cones      = n_cones
      self.n_cone_comps = n_cone_comps
      self.n_bc_comps   = n_dbp1_comps + n_hbp1_comps
      self.n_gc_comps   = n_gca_comps + n_gcaoff_comps
    
      self.comp_data = comp_data

    return comp_data
  
  ############################################################################
  def add_comsol_xyz_to_comp_data(self, verbose=True):

    # Get retsim x, y, z
    x = self.comp_data['x'].values
    y = self.comp_data['y'].values
    z = self.comp_data['z'].values
    
    if x.size == 0:
      print('Warning: Compartment data has zero length.\n'
            'Maybe you use an experiment file, that was not compiled yet.\n'
            'Fix this by (creating and) compiling the retsim expt_file.\n'
            'If this does not work your experiment might be broken.')
      return

    if verbose:
      print('\nRaw compartment data:')
      print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
      print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
      print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)))
    
    # Rotate.
    x, y, z = self.rotate_xyz(x, y, z)
    
    if verbose:
      print('\nCompartment data after rotation:')
      print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
      print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
      print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)))
    
    if self.auto_correct_comp_data:
      # Subtract comp mean.
      bp_idxs = (self.comp_data['cell_type'] == DBP1) | (self.comp_data['cell_type'] == HBP1)
      comp_x = self.comp_data['x'][bp_idxs].values
      comp_y = self.comp_data['y'][bp_idxs].values
      comp_z = self.comp_data['z'][bp_idxs].values
      comp_x, comp_y, comp_z = self.rotate_xyz(comp_x, comp_y, comp_z)

      x -= (comp_x.max()+comp_x.min())*0.5
      y -= (comp_y.max()+comp_y.min())*0.5
      z -= (comp_z.max()+comp_z.min())*0.5
      
      if verbose:
        print('\nCompartment data after subtracting BP compartment mean:')
        print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
        print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
        print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)))
      
      # Add morph mean.     
      morph_x = self.morph_data['x'].values
      morph_y = self.morph_data['y'].values
      morph_z = self.morph_data['z'].values
      morph_x, morph_y, morph_z = self.rotate_xyz(morph_x, morph_y, morph_z)
      
      x += (morph_x.max() + morph_x.min())*0.5
      y += (morph_y.max() + morph_y.min())*0.5
      z += (morph_z.max() + morph_z.min())*0.5
    
    # Add offset.
    x, y, z = self.offset_xyz(x, y, z)
    
    # Print information.
    if verbose:
      print('\nCOMSOL compartment data before scaling.')
      print('x in ({:.2f}, {:.2f})'.format(np.min(x), np.max(x)), end='\t')
      print('y in ({:.2f}, {:.2f})'.format(np.min(y), np.max(y)), end='\t')
      print('z in ({:.2f}, {:.2f})'.format(np.min(z), np.max(z)), end='\n\n')
    
    self.comp_data['comsol_x'] = x
    self.comp_data['comsol_y'] = y
    self.comp_data['comsol_z'] = z
  
  ############################################################################
  def create_comsol_comp_file(
      self, comsol_scale_factor=1e-6, verbose=True,
      inc_cones=True, inc_bc=True, inc_gc=False,
      to_file=True,
    ):  
    
    # Get indexes.
    cone_idxs = (self.comp_data['cell_type'] == XCONE)
    bp_idxs   = (self.comp_data['cell_type'] == DBP1) | (self.comp_data['cell_type'] == HBP1)
    gc_idxs   = (self.comp_data['cell_type'] == GCA) | (self.comp_data['cell_type'] == GCAOFF)
    
    # Get valid compartments.
    valid_idxs = np.zeros(self.comp_data.shape[0], dtype=bool)
    if inc_cones: valid_idxs = valid_idxs | cone_idxs
    if inc_bc:    valid_idxs = valid_idxs | bp_idxs
    if inc_gc:    valid_idxs = valid_idxs | gc_idxs
    
    # Get x, y, z
    x = self.comp_data['comsol_x'][valid_idxs].values
    y = self.comp_data['comsol_y'][valid_idxs].values
    z = self.comp_data['comsol_z'][valid_idxs].values

    # Make data frame.
    comsol_comps = pd.DataFrame(
      {'x': x*comsol_scale_factor, 'y': y*comsol_scale_factor, 'z': z*comsol_scale_factor}
    )

    if to_file:
      data_utils.make_dir('Neurons')
      path = 'Neurons/' + self.comsol_compfile
        
      # Save.
      comsol_comps.to_csv(path, index=False, header=False, columns=['x', 'y', 'z'], sep=' ')
      if verbose: print('Created compartment file @ ' + path + ' with ' + str(comsol_comps.shape[0]) + ' compartments')
    else:
      return comsol_comps
      
  ############################################################################
  def set_rot(self, mxrot=None, myrot=None, mzrot=None):
    if mxrot is not None: self.mxrot = mxrot
    if myrot is not None: self.myrot = myrot
    if mzrot is not None: self.mzrot = mzrot
      
  ############################################################################
  def set_comsol_offset(self, x0=None, y0=None, z0=None):
    if x0 is not None: self.comsol_x0 = x0
    if y0 is not None: self.comsol_y0 = y0
    if z0 is not None: self.comsol_z0 = z0
  
  ############################################################################
  def offset_xyz(self, x, y, z):
    x = np.array(x).flatten() + self.comsol_x0
    y = np.array(y).flatten() + self.comsol_y0
    z = np.array(z).flatten() + self.comsol_z0
    return x,y,z
  
  ############################################################################
  def rotate_xyz(self, x, y, z):
    
    xyz = np.vstack([np.array(x).flatten(), np.array(y).flatten(), np.array(z).flatten()])
    
    # Rotate around x.
    xyz = np.array([[1, 0, 0],
                    [0, np.cos(self.mxrot*np.pi/180), -np.sin(self.mxrot*np.pi/180)],
                    [0, np.sin(self.mxrot*np.pi/180), np.cos(self.mxrot*np.pi/180)]]) @ xyz
    
    # Rotate around y.    
    xyz = np.array([[np.cos(self.myrot*np.pi/180), 0, np.sin(self.myrot*np.pi/180)],
                    [0, 1, 0],
                    [-np.sin(self.myrot*np.pi/180), 0, np.cos(self.myrot*np.pi/180)]]) @ xyz
    
    # Rotate around z.                
    xyz = np.array([[np.cos(self.mzrot*np.pi/180), -np.sin(self.mzrot*np.pi/180), 0],
                    [np.sin(self.mzrot*np.pi/180), np.cos(self.mzrot*np.pi/180), 0],
                    [0, 0, 1]]) @ xyz
    
    return (xyz[0,:], xyz[1,:], xyz[2,:])

  ############################################################################
  def get_t_rng(self):
    return self.__t_rng
  
  ############################################################################
  def __get_rec_type(self, rec_type):
    if rec_type is not None:
      assert rec_type in self.rec_data.keys()
    elif len(self.rec_data) == 1:
      rec_type = list(self.rec_data.keys())[0]
    else:
      raise ValueError('please define rec_type')
      
    return rec_type
  
  ############################################################################
  def __get_any_Vm_V_pair(self, rec_type):
    for rec_name_Vm in self.rec_data[rec_type]['Data'].columns:
      if 'Vm' in rec_name_Vm:
        for rec_name_V in self.rec_data[rec_type]['Data'].columns:
          if rec_name_Vm.replace('Vm', 'V') in rec_name_V:
            return rec_name_Vm, rec_name_V
    return None, None
  
  ############################################################################
  def get_stim_data(self, rec_type=None):
    if self.stim_type in ['Light']:
      rec_type = self.__get_rec_type(rec_type)
      stim_data = {
        'Time': self.rec_data[rec_type]['Time'],
        'Stim': self.rec_data[rec_type]['Stim'],
        'Label': 'Light'
      }
    
    elif self.stim_type in ['Vext']:
      
      # Find a pair with Vm an V to compute Vext.
      rec_name_Vm, rec_name_V = self.__get_any_Vm_V_pair(rec_type)

      if (rec_name_Vm is not None) and (rec_name_V is not None):
        stim_Vext = self.rec_data[rec_type]['Data'][rec_name_V] - self.rec_data[rec_type]['Data'][rec_name_Vm]
        label = rec_name_Vm.replace('Vm', 'Vext')
      else:
        stim_Vext = np.full(self.rec_data[rec_type]['Time'].size, np.nan)
        warnings.warn('Have not found stimulus')
        label = 'None'
      
      stim_data = {
        'Time': self.rec_data[rec_type]['Time'],
        'Stim': stim_Vext,
        'Label': label,
      }
    
    return stim_data
  
  ############################################################################
  def get_rec_data(self, rec_type=None):
    rec_type = self.__get_rec_type(rec_type)
    return self.rec_data[rec_type]['Data']
      
  ############################################################################
  def get_rec_time(self, rec_type=None):
    if rec_type is not None:
      assert rec_type in self.rec_data.keys()
    elif len(self.rec_data) == 1:
      rec_type = list(self.rec_data.keys())[0]
    else:
      raise ValueError('please define rec_type')
      
    return self.rec_data[rec_type]['Time']
  
  ############################################################################
  def get_unit(self, param):
    if param in self.params_unit:
      return self.params_unit[param]
    else:
      return 1.0
  
  ############################################################################
  def get_bc2cone_nodes(self, set_n_cones=None, verbose=False, morphfile=None, bp_type=None, rotate=False):
    if morphfile is None:
      # If no morph file is given at least one default should be None.
      assert self.bp_morphfile_ON is None or self.bp_morphfile_OFF is None, 'Ambiguous input'
   
      # Return empty lists and use default morphology.
      if self.bp_morphfile_ON is None and self.bp_morphfile_OFF is None:
        if verbose: print('No morph file given --> No cone contacts')
        return []
      else:
        morphfile = self.bp_morphfile_ON or self.bp_morphfile_OFF
  
    compfile = 'cone_connections_' + morphfile + '.csv'
    
    bp_type = bp_type or self.bp_type
    set_n_cones = set_n_cones or self.set_n_cones
    
    # Run experiment.
    im,err_out = python2retsim.retsim(
      d=3,
      expt='get_cone_connections',
      cone_densfile='dens_cone_emtpy.n',
      dbp1_file=morphfile,
      hbp1_file=morphfile,
      nvalfile='nval_get_cone_connections.n',
      is_OFF_bp=bp_type in ['OFF', 'CBC1', 'CBC2', 'CBC3a', 'CBC3b', 'CBC4'],
      print_err=True, ninfo=3,
      mxrot=self.mxrot,
      myrot=self.myrot,
      mzrot=self.mzrot,
      cone_soma_z=self.cone_soma_z,
      set_n_cones=set_n_cones,
      timeout=self.timeout,
      print_comps=1,
      compfile='comp_files/' + compfile,
    )
    
    if verbose:
      plt.figure(figsize=(12,12))
      plt.imshow(im)
      plt.show()
    if verbose > 1: print(err_out)
    
    # Get cone contact data.
    comp_data = self.read_retsim_compfile(compfile=compfile, verbose=False, update=False)
    
    # Extract connection nodes. Assert that nodes were not erased.
    bc2cone_nodes = []
    for line in err_out.split('\n'):
      if "connect_cell prect cone" in line:
        cone = int(line[line.find("precn")+6:line.find(" prenode")])
        if cone in comp_data['cell_number'][comp_data['cell_type'] == XCONE].values:
          bc2cone_nodes.append(int(line[line.find("postnode ") + 9:line.find(" mindist ")]))
    
    assert len(bc2cone_nodes) == np.unique(comp_data['cell_number'][comp_data['cell_type'] == XCONE]).size
    if verbose and len(bc2cone_nodes) == 0: print('No cone contacts founds')
    
    return bc2cone_nodes
  
  ############################################################################
  def get_bc_node_retsim_xyz(self, nodes, morphfile=None):
  
    # Get x, y and z in retsim coordinates.
    if morphfile is None:
      morph_data = self.morph_data
    else:
      morph_data = self.read_retsim_morphfile(
        morphfile=morphfile, verbose=False, update=False
      )
    
    assert morph_data is not None
    
    x_nodes = []
    y_nodes = []
    z_nodes = []
    
    for node in nodes:
      node_idx = np.where(morph_data['node'] == node)[0]
      if node_idx.size != 1:
        print('Node not found or multiple times found: ' + str(node))
      else:
        node_idx = node_idx[0]
        x_nodes.append(morph_data['x'][node_idx])
        y_nodes.append(morph_data['y'][node_idx])
        z_nodes.append(morph_data['y'][node_idx])
  
    x_nodes = np.array(x_nodes)
    y_nodes = np.array(y_nodes)
    z_nodes = np.array(z_nodes)
  
    return x_nodes, y_nodes, z_nodes
  
  ############################################################################
  def node2comp(self, node, cell_type_idx=None):
    cell_type_idx = cell_type_idx or self.bp_type_idx
    
    for comp, nodes in self.comp_data['nodes'][self.comp_data['cell_type'] == cell_type_idx].iteritems():
      if node in nodes:
        return comp
    
    return -1
    
  ############################################################################
  def comp2node(self, comp):   
    return self.comp_data.loc[comp]['nodes'][0]
  
  ############################################################################
  def init_retsim(
    self, d=3, R=False, verbose=False, sim_params={}, print_comps=True,
    update=True, ninfo=1, expt_idx=0, **kwargs
  ):
  
    if update: assert print_comps, 'Do not update comps when you do not print' 
  
    if len(kwargs) > 0:
      print('I don\'t not what these parameters do. They will be passed directly to retsim = ', str(kwargs))
  
    retsim_params = self.get_retsim_base_params(print_comps, ninfo, expt_idx)
    retsim_params.update(self.params_default)
    retsim_params.update(kwargs)
    retsim_params.update(sim_params)
    retsim_params = self.add_units_to_params(retsim_params)
    retsim_params = self.add_adaptive_cpl_params(retsim_params)
    if verbose > 1: self.print_params(retsim_params)

    self.find_files(retsim_params)

    # Run experiment.
    im, err_out = python2retsim.retsim(d=d, R=R, **retsim_params, print_err=1)
    
    # Print error?
    if verbose: print(err_out)
    
    # Update cell image.
    self.retsim_im = im
    
    # Update comp data of cell?
    if update:
      self.read_retsim_morphfile(update=True, verbose=False)
      self.add_comsol_xyz_to_morph_data(verbose=False)
      self.read_retsim_compfile(verbose=False, update=True)
      self.add_comsol_xyz_to_comp_data(verbose=False)
    
    return im

  ############################################################################
  def get_retsim_base_params(self, print_comps, ninfo, expt_idx):
    retsim_params = {}
    retsim_params['expt'] = self.expt_file_list[expt_idx]
    
    if self.cone_densfile    is not None: retsim_params['cone_densfile']  = self.cone_densfile
    if self.bp_densfile_ON   is not None: retsim_params['dbp1_densfile']  = self.bp_densfile_ON
    if self.bp_densfile_OFF  is not None: retsim_params['hbp1_densfile']  = self.bp_densfile_OFF
    if self.bp_morphfile_ON  is not None: retsim_params['dbp1_file']      = self.bp_morphfile_ON
    if self.bp_morphfile_OFF is not None: retsim_params['hbp1_file']      = self.bp_morphfile_OFF
    if self.nval_file        is not None: retsim_params['nvalfile']       = self.nval_file
    if self.chanparams_file  is not None: retsim_params['chanparamsfile'] = self.chanparams_file

    retsim_params['print_comps'] = int(print_comps)
    retsim_params['compfile']    = os.path.join('comp_files', self.compfile)
    
    retsim_params['is_OFF_bp'] = int(self.is_OFF_bp)
    retsim_params['make_dbp1'] = int(self.make_dbp1)
    retsim_params['make_hbp1'] = int(self.make_hbp1)

    retsim_params['cone_soma_z'] = self.cone_soma_z
    retsim_params['make_cones']  = int(self.make_cones)  
    retsim_params['merge_cones'] = int(self.merge_cones)  
    retsim_params['timeout']     = self.timeout

    retsim_params['mxrot'] = self.mxrot
    retsim_params['myrot'] = self.myrot
    retsim_params['mzrot'] = self.mzrot

    retsim_params['ninfo'] = ninfo

    return retsim_params
  
  ############################################################################
  def find_files(self, retsim_params):
  
    def find_file(filename):
      if not os.path.isfile(filename):
        print('File does not exist:')
        print(filename)

    for k, v in retsim_params.items():
      if 'densfile' in k:
        if not ('cone' in k) or self.make_cones:
            find_file(os.path.join(self.retsim_path, 'runconf', v)) 
      elif 'bp1_file' in k:
        find_file(os.path.join(self.retsim_path, 'runconf', v)) 
      elif 'nvalfile' in k:
        find_file(os.path.join(self.retsim_path, 'runconf', v)) 
      elif 'chanparamsfile' in k:
        find_file(os.path.join(self.retsim_path, 'runconf', v))
  
  ############################################################################
  def add_restim_stim_params(self, retsim_params, stim_idx, rec_type, stim_type):
    assert self.syn_dt <= self.syn_dt
    assert self.sim_dt <= self.stim_dt
    assert self.sim_dt <= self.rec_dt
    assert self.sim_dt <= self.predur_dt 

    retsim_params['stimfile']  = self.create_retsim_stim_file_name(stim_idx=stim_idx)
    
    retsim_params['rectype']  = rec_type or self.rec_type
    retsim_params['stimtype'] = stim_type or self.stim_type

    retsim_params['predur']         = self.predur
    retsim_params['stimdur']        = self.stimdur
    retsim_params['set_ploti']      = self.rec_dt
    retsim_params['set_timinc']     = self.sim_dt
    retsim_params['predur_timinc']  = self.predur_dt
    retsim_params['set_stiminc']    = self.syn_dt
    retsim_params['set_srtimestep'] = self.stim_dt
    
    return retsim_params
  
  ############################################################################
  def add_units_to_params(self, retsim_params):
    for param, param_value in retsim_params.items():
      if param in self.params_unit:
        retsim_params[param] = param_value*self.get_unit(param)
    return retsim_params
    
  ############################################################################
  @staticmethod
  def print_params(retsim_params):
    print('### Run simulation with: ###')
    for k, v in retsim_params.items():
      if isinstance(v, float):  print(k.ljust(30) + " = {:.4g}".format(v))
      else:                     print(k.ljust(30) + " = " + str(v))
    print("\n")
 
   ############################################################################
  @staticmethod
  def add_adaptive_cpl_params(retsim_params):
    for param, param_value in retsim_params.items():
      if ('cpl' in param) and isinstance(param_value, dict):
        retsim_params[param] = param_value['c'] * np.sqrt(retsim_params[param_value['ri']] / retsim_params[param_value['rm']])
    return retsim_params
 
  ############################################################################
  def run(
    self, sim_params={}, rec_type=None, stim_type=None, verbose=False, reset_retsim_stim=True,
    ninfo=1, update_cell_rec_data=False, warning_level=0, DEBUG=False,
    plot=False, stim_idx=None, print_comps=False, expt_idx=0, print_params=False, **kwargs
  ):
    assert isinstance(sim_params, dict)
    
    if len(kwargs) > 0:
      print('I don\'t not what these parameters do. They will be passed directly to retsim = ', str(kwargs))
    
    if stim_idx is not None:
      assert isinstance(self.stim, list)
      assert len(self.stim) > stim_idx
    
    # Initialize retsim?
    if self.retsim_im is None: self.init_retsim()
    
    # Reset neuronc stimulus. If not check if file exists.
    if reset_retsim_stim and self.stim is not None:
      self.create_retsim_stim_file(stim_idx=stim_idx)
    else:
      err = 'Stimulus does not exist. Set reset_retsim_stim to True .'
      assert os.path.isfile(self.create_retsim_stim_file_name(stim_idx=stim_idx)), err
      
    # Set general retsim default parameters.
    retsim_params = self.get_retsim_base_params(print_comps, ninfo, expt_idx)
    retsim_params = self.add_restim_stim_params(retsim_params, stim_idx, rec_type, stim_type)
    retsim_params.update(self.params_default)
    retsim_params.update(kwargs)
    retsim_params.update(sim_params)
    retsim_params = self.add_units_to_params(retsim_params)
    retsim_params = self.add_adaptive_cpl_params(retsim_params)
    if print_params: self.print_params(retsim_params)
    
    # Run simulation.
    rec_data, labels, err_out = python2retsim.retsim(print_err=True, print_labels=True, **retsim_params)
    
    # Print Error?
    if verbose: print(err_out)
    
    # If debug, stop here.
    if DEBUG: return rec_data, labels, err_out
    
    # Check output data. If it does not pass, return Nones.
    if not self.check_retsim_output(rec_data=rec_data, labels=labels, retsim_params=retsim_params, warning_level=warning_level):
      return (None, None, None)
    
    # Rename columns of rec_data.
    if self.stim_type == "Light":
      try:
        rec_data.columns = ['Time', 'Stim'] + labels[1:]
      except:
        rec_data.columns = ['Time', 'Stim'] + [str(i) for i in range(len(rec_data.columns)-2)]
        print('label size did not match rec_data size. Labels were:' + str(labels))
    else:
      try:
        rec_data.columns = ['Time'] + labels
      except:
        rec_data.columns = ['Time'] + [str(i) for i in range(len(rec_data.columns)-1)]
        print('label size did not match rec_data size. Labels were:' + str(labels))
    
    # Get time.
    rec_time = rec_data['Time'].values
    rec_data = rec_data.drop(['Time'], axis=1)

    #if (self.syn_dt != 1e-4):
    #  rec_data.iloc[:,['rate' in col for col in rec_data.columns]] *= (self.syn_dt / 1e-4)

    # Get stimulus.
    if 'Stim' in rec_data:
      rec_stim = rec_data['Stim'].values
      rec_data = rec_data.drop('Stim', axis=1)
    else:
      rec_stim = None
    
    # Drop last point?
    t_rng = self.get_t_rng()
    
    if rec_time[-1] > (t_rng[1]-t_rng[0]):
      assert rec_time[-2] <= (t_rng[1]-t_rng[0]), str(rec_time[-2]) + ' > ' + str((t_rng[1]-t_rng[0]))
      
      rec_time = rec_time[:-1]
      if rec_stim is not None: rec_stim = rec_stim[:-1]
      rec_data = rec_data.iloc[:-1,:]
    
    # Save data.
    if update_cell_rec_data:
      self.rec_data[retsim_params['rectype']] = {} 
      self.rec_data[retsim_params['rectype']]['Stim'] = rec_stim
      self.rec_data[retsim_params['rectype']]['Time'] = rec_time
      self.rec_data[retsim_params['rectype']]['Data'] = rec_data
    
    if plot:
      plot_cell_rec_data.plot_rec_data(rec_time=rec_time, rec_data=rec_data, rec_stim=rec_stim)
    
    return (rec_data, rec_time, rec_stim)

  ############################################################################
  def check_retsim_output(self, rec_data, labels, retsim_params, warning_level=1):
    # Check if model output is valid.
    if not(isinstance(rec_data, pd.DataFrame)) or not(isinstance(labels, list)):
      if warning_level >= 2: warnings.warn(' Simulation crashed! Parameters: ' + str(retsim_params))
      if warning_level == 1: print('ERR(no dat)', end=' ')
      for file in os.listdir(self.retsim_path):
        if file.startswith("core."):
          try: os.remove(os.path.join(self.retsim_path, file))
          except: pass
      return False
    if rec_data.shape[0] <= 1:
      if warning_level >= 2: warnings.warn(' No Data recorded! Parameters: ' + str(retsim_params))
      if warning_level == 1: print('ERR(dat shape)', end=' ')
      return False
    try:
      if not(np.all(np.isfinite(rec_data.values))):
        if warning_level >= 2: warnings.warn(' Data has NaNs! Parameters: ' + str(retsim_params))
        if warning_level == 1: print('ERR(NaN)', end=' ')
        return False
    except:
      if warning_level >= 2: warnings.warn(' Data has non numbers! Parameters: ' + str(retsim_params))
      if warning_level == 1: print('ERR(Type)', end=' ')
      return False
    return True
    
  ############################################################################
  def run_parallel(self, sim_params_list=[{}], n_parallel=20, plot=True, verbose=False, legend=None):
    assert isinstance(sim_params_list, list), 'must be a list of dictionaries'
    assert len(sim_params_list) > 0, 'empty list, use empty dicts in list to use no params'
    
    if verbose:
      for idx, sim_params in enumerate(sim_params_list):
        print(str(idx+1) + ': ' + str(sim_params), end='\n\n')
    
    with Pool(processes=n_parallel) as pool:
      rec_data_list = pool.map(self.run, sim_params_list)
    
    if plot:
      plot_cell_rec_data.plot_rec_data_list(rec_data_list, sim_params_list, legend=legend)
    
    return rec_data_list
    
  ############################################################################
  def run_parallel_stimuli(self, n_parallel=20, plot=True):
    
    n_stims = len(self.stim)
    
    with Pool(processes=n_parallel) as pool:
      rec_data_list = pool.map(self.run_stim_idx, np.arange(0,n_stims))
    
    if plot:
      plot_cell_rec_data.plot_rec_data_list(
        rec_data_list, ["stim" + str(idx) for idx in range(n_stims)]
      )
    
    return rec_data_list

  ############################################################################
  def run_stim_idx(self, stim_idx=None):
    return self.run(stim_idx=stim_idx)
    
  ###########################################################################################
  def create_rd_morph(self, region='R1'):
    
    if self.bp_type_idx == DBP1: morphfile = self.bp_morphfile_ON
    else:                        morphfile = self.bp_morphfile_OFF 

    # Load morph.
    morph_data = self.read_retsim_morphfile(
      morphfile=morphfile, verbose=False, update=False
    )
    
    # Find dendritic tip nodes. Nodes that have no child but themselves.
    start_n_idx = []
    for node_idx in np.where(morph_data["region"] == 'R1')[0]:
      if node_idx in morph_data["c_idx"].iloc[node_idx]:
        start_n_idx.append(node_idx)
    
    remove_idx = []
    for n_idx in start_n_idx:
      in_region = True
      while (in_region):
        p_idx = morph_data["p_idx"][n_idx][0]
        if morph_data["region"][p_idx] == region:
          in_region = True
          remove_idx += [p_idx]
        else:
          in_region = False
            
        n_idx = p_idx
    
    remove_idx += start_n_idx
    
    print('Removing {:d} nodes'.format(len(remove_idx)))
    
    # Create new morph without dendritic nodes.
    morph_rd = morph_data[~morph_data.index.isin(remove_idx)][["node", "parent", "dia", "x", "y", "z", "region", "dend"]]

    # Save morphology.
    morph_rd.to_csv(
      self.retsim_path + 'runconf/' + morphfile + '_rd', sep='\t', index=False,
      header =["# node", "parent", "dia", "xbio", "ybio", "zbio", "region", "dendr"],
      columns=["node", "parent", "dia", "x", "y", "z", "region", "dend"]
    )
    
    # Add rd to morphfile
    if morphfile[-3:] != '_rd':
      morphfile += '_rd'
    
    # Update cell morph_file
    if self.bp_type_idx == DBP1: self.bp_morphfile_ON = morphfile
    else:                        self.bp_morphfile_OFF = morphfile 
  
  ############################################################################ 
  def save_init_data(self, dirname):
    data_utils.make_dir(dirname)
    data_utils.save_var(self.rec_data, file=dirname + '/cell_rec_data.pkl')
      
  ############################################################################ 
  def load_init_data(self, dirname):
    self.rec_data = data_utils.load_var(dirname + '/cell_rec_data.pkl')

############################################################################  
  def plot_cpl_list(self, cpl_list):
    # Compute number of compartments.
    retsim_im_list = []
    n_comp_list = []
    
    # Get update cpl and plot.
    for cpl in cpl_list:
      if isinstance(cpl, dict):
        self.update_cpl(**cpl)
      else:
        self.update_cpl(cpl)

      self.init_retsim(Verbose=False)
        
      # Extract data.
      retsim_im_list.append(self.retsim_im)
      n_comp_list.append(self.n_bc_comps)
  
    # Reset.
    cpl = cpl_list[0]
    if isinstance(cpl, dict):
      self.update_cpl(**cpl, verbose=False)
    else:
      self.update_cpl(cpl, verbose=False)
      
    # Plot.
    n_rows = int(np.ceil(len(cpl_list)/2))
    
    plt.figure(figsize=(12,n_rows*6))
    for idx, retsim_im in enumerate(retsim_im_list):
      ax = plt.subplot(n_rows, 2, idx+1)
      plt.imshow(retsim_im)
      plt.axis('off')
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      plt.title(str(n_comp_list[idx]))

############################################################################
############################################################################

class Cone(Cell):
  
  def __init__(
      self, stimulus=None, t_rng=None, timeout=60*60*1,
      retsim_path='/gpfs01/berens/user/joesterle/berens/neuronc-master/nc/models/retsim/',
      **kwargs
    ):
    super().__init__(timeout=timeout, retsim_path=retsim_path)
    
    self.stim_type = 'Light'
    self.bp_morphfile_ON = 'morph_single_comp'
    
    self.read_retsim_morphfile()
    
    self.cone_densfile = 'dens_empty.n'
    self.nval_file = 'nval_default.n'
    self.chanparams_file = 'chanparams_empty.n'
    self.expt_base_file_list = ['retsim_files/expt_optimize_cones_base.cc']
    self.expt_file_list = ['optimize_cones']
    self.retsim_stim_file_base = 'Light_stimulus_optimize_cones'
    
    self.compfile = 'comp_single_comp.csv'
    self.make_cones = True
    self.set_n_cones = None
    self.cone_soma_z = 20
    self.bp_type_idx = DBP1
    self.gc_type_idx = -1
    self.is_OFF_bp = False
    self.make_dbp1 = True
    self.make_hbp1 = False
    
    # Initialize other values.
    self.__dict__.update(kwargs)
    
    if stimulus is not None: self.set_stim(stimulus, self.stim_type)
    if stimulus is not None or t_rng is not None: self.update_t_rng(t_rng)
  
############################################################################
class CBC(Cell):
  
  def __init__(
    self, bp_type,
    make_cones=True, merge_cones=True,
    stim_type='Light', stimulus=None,
    t_rng=None, timeout=60*60*5, predur=1,
    expt_file=None, expt_file_list=None,
    expt_base_file=None, expt_base_file_list=None,
    bp_densfile=None, bp_morphfile=None,
    cone_densfile=None,
    nval_file=None, chanparams_file=None,
    compfile=None, comsol_compfile=None,
    retsim_stim_file_base=None,
    params_default={}, params_unit={},
    retsim_path='/gpfs01/berens/user/joesterle/berens/neuronc-master/nc/models/retsim/',
  ):
    if expt_file is not None:
      assert expt_file_list is None
      expt_file_list = [expt_file]     
      
    if expt_base_file is not None:
      assert expt_base_file_list is None
      expt_base_file_list = [expt_base_file]
  
    super().__init__(timeout=timeout, retsim_path=retsim_path)
    
    # Set type to ON or OFF.
    if bp_type in ['CBC1', 'CBC2', 'CBC3a', 'CBC3b', 'CBC4']:
      self.bp_type_idx = HBP1
      self.gc_type_idx = GCAOFF
      self.myrot = 20
      self.is_OFF_bp = True
      self.make_dbp1 = False
      self.make_hbp1 = True
    elif bp_type in ['CBC5t', 'CBC5o', 'CBC5i', 'CBCX', 'CBC6', 'CBC7', 'CBC8', 'CBC9']:
      self.bp_type_idx = DBP1
      self.gc_type_idx = GCA
      self.myrot = 90
      self.is_OFF_bp = False
      self.make_dbp1 = True
      self.make_hbp1 = False
    else:
      raise ValueError('Not Implemented')
    
    self.bp_type = bp_type
    
    # Set number of cones.
    if   bp_type == 'CBC1':  self.set_n_cones = 4
    elif bp_type == 'CBC2':  self.set_n_cones = 4
    elif bp_type == 'CBC3a': self.set_n_cones = 5
    elif bp_type == 'CBC3b': self.set_n_cones = 5
    elif bp_type == 'CBC4':  self.set_n_cones = 6
    elif bp_type == 'CBC5t': self.set_n_cones = 3
    elif bp_type == 'CBC5o': self.set_n_cones = 3
    elif bp_type == 'CBC5i': self.set_n_cones = 7
    elif bp_type == 'CBCX':  self.set_n_cones = 2
    elif bp_type == 'CBC6':  self.set_n_cones = 4
    elif bp_type == 'CBC7':  self.set_n_cones = 5
    elif bp_type == 'CBC8':  self.set_n_cones = 7
    elif bp_type == 'CBC9':  self.set_n_cones = 4
    else: ValueError('Not Implemented')
  
    # Set morph_file.
    if   bp_type == 'CBC1':  self.bp_morphfile_OFF = bp_morphfile or 'cbp_0395_t1'
    elif bp_type == 'CBC2':  self.bp_morphfile_OFF = bp_morphfile or 'cbp_0395_t2'
    elif bp_type == 'CBC3a': self.bp_morphfile_OFF = bp_morphfile or 'cbp_0476_t3b'
    elif bp_type == 'CBC3b': self.bp_morphfile_OFF = bp_morphfile or 'cbp_0476_t3b'
    elif bp_type == 'CBC4':  self.bp_morphfile_OFF = bp_morphfile or 'cbp_0510_t4'
    elif bp_type == 'CBC5t': self.bp_morphfile_ON  = bp_morphfile or 'cbp_0572_t5t'
    elif bp_type == 'CBC5o': self.bp_morphfile_ON  = bp_morphfile or 'cbp_0576_t5o'
    elif bp_type == 'CBC5i': self.bp_morphfile_ON  = bp_morphfile or 'cbp_0540_t5i'
    elif bp_type == 'CBCX':  self.bp_morphfile_ON  = bp_morphfile or 'cbp_0592_t5x'
    elif bp_type == 'CBC6':  self.bp_morphfile_ON  = bp_morphfile or 'cbp_0647_t6'
    elif bp_type == 'CBC7':  self.bp_morphfile_ON  = bp_morphfile or 'cbp_0658_t7'
    elif bp_type == 'CBC8':  self.bp_morphfile_ON  = bp_morphfile or 'cbp_0686_t8'
    else: ValueError('Not Implemented')
  
    # Set channel density file.
    if self.is_OFF_bp:
      self.bp_densfile_OFF = bp_densfile or 'dens_empty.n'
    else:
      self.bp_densfile_ON = bp_densfile or 'dens_empty.n'
    
    self.cone_densfile = cone_densfile or 'dens_empty.n'
    self.nval_file = nval_file or 'nval_default.n'
    
    # Set experiment file(s).
    self.expt_base_file_list = expt_base_file_list or ['retsim_files/expt_optimize_bc_base.cc'] # Base file.
    self.expt_file_list      = expt_file_list      or [bp_type + '_optimize'] # Output file.
    self.stim_type = stim_type
    self.predur = predur
    self.make_cones = make_cones
    self.merge_cones = merge_cones

    # Output files.
    self.retsim_stim_file_base = retsim_stim_file_base or f'{stim_type}_{bp_type}_optimize'
    
    # Set parameters.
    self.params_default = params_default.copy()
    self.params_unit    = params_unit.copy()

    # Set chanparams file.
    self.compfile        = compfile        or 'comp_'        + bp_type + '_' + str(np.random.randint(0,100)) + '.csv'
    self.comsol_compfile = comsol_compfile or 'comsol_comp_' + bp_type + '_' + str(np.random.randint(0,100)) + '.txt'
    self.chanparams_file = chanparams_file or 'chanparams_empty.n'
    
    self.read_retsim_morphfile()
    
    # Update cpl.
    try:
      self.update_cpl(**params_default, verbose=False)
    except:
      print('Warning: Could not initialize cell.\n'
            'Maybe you use an experiment file, that was not compiled yet.\n'
            'Fix this by (creating and) compiling the retsim expt_file.\n'
            'If this does not work your experiment might be broken.')
    
    # Set stimulus and t_rng.
    if stimulus is not None: self.set_stim(stimulus, stim_type)
    if stimulus is not None or t_rng is not None: self.update_t_rng(t_rng)
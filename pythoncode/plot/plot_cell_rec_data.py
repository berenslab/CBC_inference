import numpy as np
from matplotlib import pyplot as plt
import warnings

############################################################################
def plot_rec_data(rec_time, rec_data, rec_stim=None):
  
  ''' Plot recorded data.
  
  Parameters:
  
  rec_time : 1d-array
    Time of recording.
    
  rec_data : DataFrame
    All recorded parameters in columns.
    Number of rows should match rec_time size.
    
  rec_stim : 1d-array or None
    Recorded stimulus.
  '''
  
  # All column names.
  data_names = list(np.unique(rec_data.columns))
  
  # Get list of cases, where Vext can be computed from V and Vm.
  Vext_refs = get_Vext_refs(data_names)
  
  # Plot.
  n_subplots = (rec_stim is not None) + len(Vext_refs) + len(data_names)
  fig, axs = plt.subplots(n_subplots, 1, figsize=(12,1.4*n_subplots), sharex=True, squeeze=False)
  
  axs = axs.flatten()
  
  ax_i = 0
  
  # Plot stimulus.
  if rec_stim is not None:    
    axs[ax_i].plot(rec_time, rec_stim)
    axs[ax_i].set_title('Stimulus')
    axs[ax_i].set_ylabel('[P*]')
    ax_i += 1
  
  # Plot Vext, if possible.
  ax_is = np.arange(ax_i, ax_i+len(Vext_refs))
  for ax_i, Vext_ref in zip(ax_is, Vext_refs):
    axs[ax_i].plot(rec_time, 1e3*(rec_data[Vext_ref['V']].values - rec_data[Vext_ref['Vm']].values))
    axs[ax_i].set_title(Vext_ref['Vext'])
    axs[ax_i].set_ylabel('[mV]')
  
  if len(Vext_refs) > 0: ax_i += 1
  
  # Plot others.
  ax_is = np.arange(ax_i, ax_i+len(data_names))
  for ax_i, data_name in zip(ax_is, data_names):
    axs[ax_i].set_title(data_name)
    if "V" in data_name:
      axs[ax_i].plot(rec_time, 1e3*rec_data[data_name])
      axs[ax_i].set_ylabel('[mV]')
    elif "Ca" in data_name:
      axs[ax_i].plot(rec_time, 1e6*rec_data[data_name])
      axs[ax_i].set_ylabel('[mM]')
    else:
      axs[ax_i].plot(rec_time, rec_data[data_name])
      axs[ax_i].set_ylabel('[?]')
  
  # Show.
  axs[-1].set_xlabel('Time [s]')
  plt.tight_layout()
  plt.show()
  
############################################################################
def plot_rec_data_list(rec_data_list, sim_params_list, legend=None):
  
  ''' Plot recorded data.
    
  Parameters:
  
  rec_data_list : list of tuples with (rec_data, rec_time, rec_stim)
    Data output of cell solver.
    
  sim_params_list : list of dicts
    Simulation parameter dict corresponding to rec_data_list.
    I.e. rec_data_list[i] was simulated with sim_params_list[i]
    
  legend : list of str
    Parameters to display in legend.
    E.g. a parameters that varied between samples.
  '''
  
  cm = plt.get_cmap('gist_rainbow')
  colors = [cm(1.*rec_data_i/len(rec_data_list)) for rec_data_i in range(len(rec_data_list))] 
   
  # All column names.
  data_names = list(np.unique(rec_data_list[0][0].columns))
  
  # Get list of cases, where Vext can be computed from V and Vm.
  Vext_refs = get_Vext_refs(data_names)
  
  # Plot.
  n_subplots = len(Vext_refs) + len(data_names)
  fig, axs = plt.subplots(n_subplots, 1, figsize=(12,1.4*n_subplots), sharex=True)
  
  ax_i = 0
  
  # Plot Vext.
  ax_is = np.arange(ax_i, ax_i+len(Vext_refs))
  for ax_i, Vext_ref in zip(ax_is, Vext_refs):
    for (rec_data, rec_time, _), sim_params in zip(rec_data_list, sim_params_list):
      if rec_data is not None:
        label = get_label(legend, sim_params)
        rec_data_Vext = (rec_data[Vext_ref['V']].values - rec_data[Vext_ref['Vm']].values)
        axs[ax_i].plot(rec_time, 1e3*rec_data_Vext, label=label)
      else:
        warnings.warn('No rec_data.')
    axs[ax_i].set_title(Vext_ref['Vext'])
    axs[ax_i].set_ylabel('[mV]')
  
  if len(Vext_refs) > 0: ax_i += 1
  
  # Plot others.
  ax_is = np.arange(ax_i, ax_i+len(data_names))
  for ax_i, data_name in zip(ax_is, data_names):
    axs[ax_i].set_title(data_name)

    if "V" in data_name:
      axs[ax_i].set_ylabel('[mV]')
    elif "Ca" in data_name:
      axs[ax_i].set_ylabel('[mM]')
    else:
      axs[ax_i].set_ylabel('[?]')
    
    for (rec_data, rec_time, _), sim_params in zip(rec_data_list, sim_params_list):
    
      label = get_label(legend, sim_params)
    
      if rec_data is not None:
        if "V" in data_name:
          axs[ax_i].plot(rec_time, 1e3*rec_data[data_name], label=label)
        elif "Ca" in data_name:
          axs[ax_i].plot(rec_time, 1e6*rec_data[data_name], label=label)
        else:
          axs[ax_i].plot(rec_time, rec_data[data_name], label=label)

      else:
        warnings.warn('No rec_data.')

    if legend is not None:
      fontsize = np.clip(14-len(sim_params_list), 6, None)
      axs[0].legend(fontsize=fontsize, bbox_to_anchor=(1, 1), loc='upper left')
   
  axs[-1].set_xlabel('Time [s]')
  plt.tight_layout()
  plt.show()
  
############################################################################
def get_Vext_refs(data_names):
  ''' Compute list of cases, where Vext can be computed from V and Vm.
  
  Parameters:
  
  data_names : list of str
    Data column names.
    
  Returns:
  
  Vext_refs : list of dicts
    Elements, where Vm and V was found, and Vext can therefore be derived.
  '''

  Vext_refs = []   
  for name in data_names:
    if "Vm" in name:
      Vm_name = name
      V_name = Vm_name.replace("Vm", "V")
      if V_name in data_names:
        Vext_name = Vm_name.replace("Vm", "Vext")
        Vext_refs.append({'Vm': Vm_name, 'V': V_name, 'Vext': Vext_name})
        
  return Vext_refs
  
############################################################################
def get_label(legend, sim_params):   
  ''' Get label for sim_params and legend setting.
  '''

  if isinstance(legend, list):
    reduced_dict = {}
    for key in legend:
      if key in sim_params:
        reduced_dict[key] = sim_params[key]
      else:
        print('Key ' + str(key) + ' not in params')
    label = str(reduced_dict)
  else:
    label = str(sim_params)
  
  return label
import data_utils
from matplotlib import pyplot as plt
import math_utils
import interpolation_utils
import numpy as np
from copy import deepcopy

############################################################################
def test_cones(cell, filename, t_rng=None):

  rec_data, rec_time, rec_stim, original_model_output =\
    __run_cell(cell, filename, t_rng, use_single_bc_comp=True)
  
  # Plot comparison.
  if rec_data is not None:
    __plot_data(
      rec_data, rec_time, rec_stim, original_model_output,
      rate_name='rate Cone', Vm_name='Cone Vm Soma1'
    )

    __numerical_comparison(
      time1=original_model_output['Time'], trace1=original_model_output['rate'],
      time2=rec_time, trace2=rec_data['rate Cone']
    )

  else:
    print('rec_data was None')
    
  return rec_data, rec_time, rec_stim, original_model_output
    
############################################################################
def test_CBC(cell, filename, t_rng=None):

  rec_data, rec_time, rec_stim, original_model_output =\
    __run_cell(cell, filename, t_rng, use_single_bc_comp=False)
  
  # Plot comparison.
  if rec_data is not None:
    __plot_data(
      rec_data, rec_time, rec_stim, original_model_output,
      rate_name='rate BC', Vm_name='BC Vm Soma'
    )
    
    __numerical_comparison(
      time1=original_model_output['Time'], trace1=original_model_output['rate'],
      time2=rec_time, trace2=rec_data['rate BC'].mean(axis=1),
    )

  else:
    print('rec_data was None')

############################################################################
def __run_cell(cell, filename, t_rng, use_single_bc_comp=False):
  
  cell = deepcopy(cell)
  
  original_model_output = data_utils.load_var(filename)
  
  # Run with original stimulus.
  cell.set_stim(original_model_output['Stimulus'], stim_type='Light')
  if t_rng is None: t_rng = original_model_output['t_rng']
  cell.update_t_rng(t_rng)
  if use_single_bc_comp: cell.update_cpl(cpl=2, verbose=False)
  cell.predur = original_model_output['predur']
  print('Running with n_bc_comps =', cell.n_bc_comps)
  
  im = cell.init_retsim(verbose=False)
  
  plt.figure(figsize=(8,8))
  plt.imshow(im)
  plt.show()
  
  try:
    rec_data, rec_time, rec_stim = cell.run(rec_type='test', plot=False, verbose=True, reset_retsim_stim=True)
  except KeyboardInterrupt:
    print("KeyboardInterrupt")
    rec_data, rec_time, rec_stim = None, None, None
  except Exception as e:
    print("Error in Simulation\n", e)
    rec_data, rec_time, rec_stim = None, None, None
  
  if rec_time is not None: rec_time += t_rng[0]
  
  return rec_data, rec_time, rec_stim, original_model_output
  
############################################################################
def __plot_data(rec_data, rec_time_plot, rec_stim, original_model_output, rate_name, Vm_name):
  fig, axs = plt.subplots(3, 1, figsize=(12,4), sharex=True)
  
  axs[0].plot(rec_time_plot, math_utils.normalize(rec_stim), label='This model')
  axs[0].plot(
    rec_time_plot,
    interpolation_utils.in_ex_polate(x_old=original_model_output['Stimulus']['Time'], y_old=original_model_output['Stimulus']['Stim'], x_new=rec_time_plot),
    '--', label='Target'
  )
  
  rate_data = rec_data[rate_name]
  if rate_data.values.ndim > 1: rate_data = rate_data.mean(axis=1)
  
  axs[1].plot(rec_time_plot, rate_data)
  axs[1].plot(original_model_output['Time'], original_model_output['rate'], '--')
  
  Vm_data = rec_data[Vm_name]
  if Vm_data.values.ndim > 1: Vm_data = Vm_data.mean(axis=1)
  
  axs[2].plot(rec_time_plot, Vm_data)
  axs[2].plot(original_model_output['Time'], original_model_output['Vm'], '--')
  
  axs[0].legend()
  plt.tight_layout()
  plt.show()
  
############################################################################
def __numerical_comparison(time1, trace1, time2, trace2):
  
  time1 = np.asarray(time1)
  time2 = np.asarray(time2)
  trace1 = np.asarray(trace1)
  trace2 = np.asarray(trace2)
  
  if (time1.size == time2.size) and np.allclose(time1, time2):
    print('Times are the same.')
  
  else:
    print('Times are not the same. Use interpolation')
    time_interpol = time1[(time1>=np.max([time1[0], time2[0]])) & (time1<=np.min([time1[-1], time2[-1]]))]
  
    trace1 = interpolation_utils.in_ex_polate(x_old=time1, y_old=trace1, x_new=time_interpol)
    trace2 = interpolation_utils.in_ex_polate(x_old=time2, y_old=trace2, x_new=time_interpol)
    
  if np.all(trace1 == trace2):
    print('Traces are exactly the same.') 
  elif np.allclose(trace1, trace2, rtol=1e-2, atol=1e-2):
    max_err = np.max(np.abs(trace2-trace1))
    print('Traces are very close, differences might be due to rounding errors. Max error = {:.2g}'.format(max_err))
  else:
    max_err = np.max(np.abs(trace2-trace1))
    print('Traces are not the equal. Max error = {:.2g}'.format(max_err))
    
#############################################################################
#def test_CBCs(self, filename_ON, filename_OFF, t_rng=None):
#  ON_model_output = data_utils.load_var(filename_ON)
#  OFF_model_output = data_utils.load_var(filename_OFF)
#          
#  backup_stim = self.stim
#  backup_stim_type = self.stim_type
#  
#  # Run with original stimulus.
#  assert np.all(ON_model_output['Stimulus'] == OFF_model_output['Stimulus'])
#  self.set_stim(ON_model_output['Stimulus'], stim_type='Light')
#  
#  if t_rng is not None: self.update_t_rng(t_rng)
#  else:                 self.update_t_rng(cone_model_output['t_rng'])
#  
#  try:
#    rec_data = self.run(
#      sim_params={'g_ac_hb': 0.0, 'g_db_ac': 0.0},
#      rec_type='test', plot=False, verbose=False, reset_retsim_stim=True,
#    )
#  except:
#    rec_data = None
#    pass
#  
#  # Reset stimulus.
#  self.set_stim(backup_stim, backup_stim_type)
#  
#  # Plot comparison.
#  if rec_data is not None:
#    plt.figure(figsize=(12,6))
#    plt.subplot(511)
#    plt.plot(rec_data[1]+self.get_t_rng()[0], math_utils.normalize(rec_data[2]), label='This model')
#    plt.plot(ON_model_output['Stimulus']['Time'], ON_model_output['Stimulus']['Stim'], '--', label='Target')
#    plt.legend()
#    
#    plt.subplot(512)
#    plt.plot(rec_data[1]+self.get_t_rng()[0], rec_data[0]['rate BC ON'].mean(axis=1))
#    plt.plot(ON_model_output['Time']+ON_model_output['t_rng'][0], ON_model_output['rate'], '--')
#    
#    plt.subplot(513)
#    plt.plot(rec_data[1]+self.get_t_rng()[0], rec_data[0]['BC Vm Soma ON'])
#    plt.plot(ON_model_output['Time']+ON_model_output['t_rng'][0], ON_model_output['Vm'], '--')
#    
#    plt.subplot(514)
#    plt.plot(rec_data[1]+self.get_t_rng()[0], rec_data[0]['rate BC OFF'].mean(axis=1))
#    plt.plot(OFF_model_output['Time']+OFF_model_output['t_rng'][0], OFF_model_output['rate'], '--')
#    
#    plt.subplot(515)
#    plt.plot(rec_data[1]+self.get_t_rng()[0], rec_data[0]['BC Vm Soma OFF'])
#    plt.plot(OFF_model_output['Time']+OFF_model_output['t_rng'][0], OFF_model_output['Vm'], '--')
#    
#    plt.tight_layout()
#    plt.show()
#  else:
#    print('rec_data was None')
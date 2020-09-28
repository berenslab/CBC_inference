from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

############################################################################    
def plot_best_samples(samples, time, loss=None, n=1, figsize=(12,6), name_rate='rate', name_Vm='Vm'):
  assert name_rate in samples
  assert name_Vm in samples
  
  d_sort_index = np.argsort(samples['loss']['total'])
  
  # Plot best neuronC data.
  plt.figure(figsize=figsize)
  
  ax = plt.subplot2grid((4,1), (0,0), rowspan=1)
  plt.title('Release Rate')
  plt.plot(time, samples[name_rate][d_sort_index[0:n],].T)
  plt.yticks([samples[name_rate][d_sort_index[0:n],].min(axis=None),\
              samples[name_rate][d_sort_index[0:n],].max(axis=None)])
  plt.ylabel('ves/s/syn')
  ax.set_xticklabels([])
  
  ax = plt.subplot2grid((4,1), (1,0), rowspan=1)
  plt.title('Membrane voltage')
  plt.plot(time, samples[name_Vm][d_sort_index[0:n],].T)
  plt.yticks([samples[name_Vm][d_sort_index[0:n],].min(axis=None),\
              samples[name_Vm][d_sort_index[0:n],].max(axis=None)])
  plt.ylabel(name_Vm)
  ax.set_xticklabels([])
  
  ax = plt.subplot2grid((4,1), (2,0), rowspan=2)
  if loss is not None:
    plt.title('iGluSnFR')
    plt.plot(loss.target_time, loss.target, c='k', zorder=3, label='target', alpha=0.5)
    for idx in range(n):
      plt.plot(loss.target_time, loss.rate2best_iGluSnFR_trace(samples[name_rate][d_sort_index[0:n],][idx,])[0])
    plt.legend(loc='upper right')
  else:
    plt.title('iGluSnFR - No available without loss')
  plt.ylabel('dF/F')
  ax.set_yticks([])

  plt.xlabel('Time [s]')
  
  plt.tight_layout()
  plt.show()

############################################################################    
def plot_samples(samples, time, n_max=100, name_rate='rate'):
  if samples is None: print('No samples given!'); return
  n_samples = samples[name_rate].shape[0]
  # Plot.
  plt.figure(figsize=(12,5))
  plt.subplot(311)
  plt.ylabel('Vesicles/s')
  plt.xlabel('Time [s]')
  if n_samples > n_max: print('Plot ' + str(n_max) + ' random samples!')
  idx = np.random.choice(np.arange(0,n_samples), size=np.min([n_max, n_samples]), replace=False)
  plt.plot(time, samples[name_rate][idx,:].T)
  ax = plt.subplot(312)
  plt.ylabel('Loss')
  for loss_name, loss in samples['loss'].items():
    if loss_name != 'total':
      plt.plot(loss, '.-', label=loss_name, linewidth=0.3, markersize=1)
  plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
  ax.set_xticklabels([])
  plt.subplot(313)
  plt.ylabel('Loss')
  plt.plot(samples['loss']['total'], '.-', label='total', linewidth=0.3, markersize=1)
  plt.xlabel('# samples')
  plt.ylabel('Loss')
  plt.legend(bbox_to_anchor=(1, 1))
  plt.tight_layout()
  plt.show()
  
############################################################################
def plot_execution_time(samples, lines=[], solid_lines=[], timeout=None, verbose=True, figsize=(10,5)):
  # Analyze cpu time of samples.
  if verbose:
    print('max    (time / sample) [min] = {:.2g}'.format(np.max(samples['wall-time']) / 60))
    print('median (time / sample) [min] = {:.2g}'.format(np.median(samples['wall-time']) / 60))
  
  walltime_idx = np.argsort(samples['wall-time'])
  
  plt.figure(figsize=figsize)
  
  ax = plt.subplot(311)
  plt.title('Wall-time')
  plt.plot(np.arange(1,1+samples['wall-time'].size), samples['wall-time'], '.-')
  plt.ylabel('Time [s]')
  
  for line in lines:
    plt.axvline(line, c='r', linestyle='--')
  for line in solid_lines:
    plt.axvline(line, c='r', linestyle='-')
  
  ax = plt.subplot(312)
  plt.title('Wall-time sorted')
  plt.plot(np.arange(1,1+samples['wall-time'].size), samples['wall-time'][walltime_idx], '.-')
  plt.ylabel('Time [s]')
  
  ax = plt.subplot(313)
  plt.title('loss(Wall-Time)')
  plt.plot(samples['wall-time'][walltime_idx], samples['loss']['total'][walltime_idx], '.')
  if timeout is not None: plt.axvline(timeout, c='r')
  plt.ylabel('loss iGluSnFR')
  plt.xlabel('Time [min]')
  
  plt.tight_layout()
  plt.show()


def plot_loss_rounds(samples, n_samples, equal_x=False):
  loss_keys = samples['loss'].keys()
  fig, axs = plt.subplots(len(loss_keys), 2, figsize=(12,2*len(loss_keys)), sharex=False, sharey='row')

  for loss_key, axs_row in zip(loss_keys, axs):
    axs_row[0].set_ylabel(loss_key, rotation=0, ha='right', va='center')
    
    for i, (lbi, ubi) in enumerate(zip(np.append(0, n_samples[:-1]).astype(int), n_samples)):
      loss_data_i = samples['loss'][loss_key][lbi:ubi]
      finite_idx = np.isfinite(loss_data_i)
      
      if equal_x:
        x_data = np.linspace(0,1,loss_data_i.size)
      else:
        x_data = np.arange(0,loss_data_i.size)
      
      axs_row[0].plot(x_data[:np.sum(finite_idx)], np.sort(loss_data_i[finite_idx]), label=str(i), clip_on=False, lw=2, c='C'+str(i))
      
      if np.sum(~finite_idx) > 0:
        axs_row[0].fill_between(
          x_data[np.sum(finite_idx):],
          np.full(np.sum(~finite_idx), np.nanmin(loss_data_i)),
          np.full(np.sum(~finite_idx), np.nanmax(loss_data_i)),
          label='_', alpha=0.3, color='C'+str(i)
        )
      sns.distplot(loss_data_i, vertical=True, hist=True, ax=axs_row[1], color='C'+str(i), kde=False)
      axs_row[1].axhline(np.mean(loss_data_i), color='C'+str(i), lw=n_samples.size-i/2, alpha=1)
    
    axs_row[0].axhline(0, c='k', ls='--', alpha=0.4)
    axs_row[1].axhline(0, c='k', ls='--', alpha=0.4)
    
    if equal_x:
      axs_row[0].set_xlim(0, 1)
    else:
      axs_row[0].set_xlim(0, )
    
  axs[0,0].legend(bbox_to_anchor=(1,1), loc='upper left', fontsize=8)
    
  plt.tight_layout()
  plt.show()
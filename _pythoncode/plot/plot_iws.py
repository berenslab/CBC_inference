from matplotlib import pyplot as plt
import numpy as np
  
def plot_iws(tds, pseudo_obs_dim, pseudo_obs=None, kernel_bandwidths=None):
  fig, axs = plt.subplots(1,4,figsize=(12,3))

  for i, td_i in enumerate(tds):
    axs[0].plot(td_i[1][:,pseudo_obs_dim], td_i[2], '.', label=str(i), color='C'+str(i))

    if (pseudo_obs is not None) and (kernel_bandwidths is not None):
      axs[0].axvline(pseudo_obs[i][0,pseudo_obs_dim]+kernel_bandwidths[i][pseudo_obs_dim], color='C'+str(i))
    
  axs[0].set_xlabel('MSE loss')
  axs[0].set_ylabel('iws')  
  axs[0].set_ylim(0, None)
  axs[0].legend()

  for i, td_i in enumerate(tds):
    axs[1].plot(np.sort(td_i[2]))
    axs[1].set_xlabel('iws')

  for i, td_i in enumerate(tds):
    axs[2].plot(np.sort(td_i[1][:,pseudo_obs_dim]))
    axs[2].set_xlabel('MSE loss')
    
  for i, td_i in enumerate(tds):
    axs[3].plot(np.sort(np.sum(np.abs(td_i[1]), axis=1)))
    axs[3].set_xlabel('total loss')


  plt.tight_layout()
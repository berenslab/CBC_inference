from matplotlib import pyplot as plt
import numpy as np

def plot(obs=None, bw=None, figsize=(10,2)):
  obs = np.asarray(obs).flatten()
  bw = np.asarray(bw).flatten()
  # Plot observed.
  plt.figure(figsize=figsize)
  plt.subplot(1,2,1)
  plt.title('obs')
  if obs is not None:
    plt.plot(np.arange(1,1+obs.size), obs, '*-')
    plt.xticks(np.arange(1,1+obs.size))
  plt.xlabel('round')
  plt.subplot(1,2,2)
  plt.title('bw')
  if bw is not None:
    plt.plot(np.arange(1,1+bw.size), bw, '*-')
    plt.xticks(np.arange(1,1+bw.size))
  plt.xlabel('round')
  plt.show()
  
def plot_logs(logs, figsize=(15,1)):
  # Plot loss of all training rounds.
  plt.figure(figsize=(figsize[0],len(logs)+figsize[1]))
  for i in range(len(logs)):
      plt.subplot(len(logs),1,i+1)
      plt.plot(logs[i]['loss'], '.-')
      plt.ylabel('NN Loss')
  
  plt.tight_layout()
  plt.xlabel('Iteration')
  plt.show()
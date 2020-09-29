from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

import interpolation_utils

############################################################################
############################################################################
class Loss():
  target   = None
  rec_time = None
  t_drop   = 0
  
  ############################################################################
  def __init__(self, target, t_drop):
    self.target = target
    self.t_drop = t_drop
    self.max_loss = {'total': 4, 'single': 1}
    
  ############################################################################
  def calc_loss(self, rec_data, plot=False):
  
    loss = {}
  
    for f0, rec_data_f in rec_data.items():
      for V0, rec_data_f_i in rec_data_f.items():
        
        if rec_data_f_i is None:
          loss['f0='+str(f0)+'_V0'+str(V0)] =  self.max_loss['single']
        
        else:
        
          trace_time = rec_data_f_i['Time']
          trace      = rec_data_f_i['Current']
      
          target_time = self.target[f0][V0]['Time']
          target      = self.target[f0][V0]['Current']
          
          loss_time = np.linspace(self.t_drop, target_time.max(), 1000)
          
          trace  = interpolation_utils.in_ex_polate(x_old = trace_time,  y_old = trace,  x_new = loss_time)
          target = interpolation_utils.in_ex_polate(x_old = target_time, y_old = target, x_new = loss_time)
          
          trace_loss  = np.mean(((trace - target)/target.max())**2)
          
          if trace_loss >= self.max_loss['single']:
            trace_loss = self.max_loss['single']
          
          loss['f0='+str(f0)+'_V0'+str(V0)] = trace_loss
          
          if plot:
            plt.figure(figsize=(12,1))
            plt.title('f0='+str(f0)+'_'+str(V0)+" --> loss = {:.4g}".format(trace_loss))
            plt.plot(loss_time, trace,  label='trace')
            plt.plot(loss_time, target, label='target')
            plt.legend()
            plt.show()
    
    loss['total'] = np.sum(list(loss.values()))
    
    return loss
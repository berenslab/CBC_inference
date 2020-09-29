import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import warnings

from loss_funcs import Loss

############################################################################
class LossOptimizeStimulation(Loss):
  base_rates = None
  loss_params = {}
  
  ############################################################################
  def __init__(self, init_rec_data, cell_target, mode='factor', p=1, maxrel=None):
    self.cells = list(init_rec_data.keys())

    self.reldur = init_rec_data[self.cells[0]]['Time'][-1] - init_rec_data[self.cells[0]]['Time'][0]
    assert self.reldur == init_rec_data[self.cells[1]]['Time'][-1] - init_rec_data[self.cells[1]]['Time'][0]

    self.set_cell_target(cell_target)
    assert self.cell_target in init_rec_data.keys()

    self.set_base(init_rec_data)
    self.n_cells = len(self.base_rates)

    assert self.n_cells == 2

    self.mode = mode
    self.p = p
    
    if 'maxrel' in mode:
      assert maxrel is not None
      self.maxrel = maxrel
    
    print('Mode is: ' + self.mode)
    print('Power is: ' + str(self.p))
    
  ############################################################################
  def set_base(self, init_rec_data):
    self.base_rates = {cell: rec_data['Data']['rate'].mean() for cell, rec_data in init_rec_data.items()}
    self.base_rel = {cell: base_rate*self.reldur for cell, base_rate in self.base_rates.items()}
    
  ############################################################################
  def set_cell_target(self, cell_target):
    if isinstance(cell_target, int):
      self.cell_target = self.cells[cell_target]
    elif cell_target in self.cells:
      self.cell_target = cell_target
    else:
      raise NotImplementedError(cell_target)
    
    print('Cell target is: ' + self.cell_target)
   
  ############################################################################
  def get_rates(self, rec_data_dict):
    rates = {cell: rec_data_dict[cell]['rate'].mean() for cell in self.cells}
    return rates
    
  ############################################################################
  def calc_loss(self, rec_data_dict, verbose=False):

    rates = self.get_rates(rec_data_dict)

    if self.mode == 'factor':
      # Compute rates relative to base rate.
      factor0 = (rates[self.cells[0]] / self.base_rates[self.cells[0]])
      factor1 = (rates[self.cells[1]] / self.base_rates[self.cells[1]])

    elif self.mode == 'abs':
      # Compute absolute rates.
      factor0 = rates[self.cells[0]]
      factor1 = rates[self.cells[1]]

    elif self.mode == 'relabs':
      # Compute absolute rates.
      factor0 = rates[self.cells[0]] - self.base_rates[self.cells[0]]
      factor1 = rates[self.cells[1]] - self.base_rates[self.cells[1]]
      
    elif self.mode == 'maxrel':
      # Compute release relative to max possible release.
      factor0 = (rates[self.cells[0]]*self.reldur / self.maxrel[self.cells[0]])
      factor1 = (rates[self.cells[1]]*self.reldur / self.maxrel[self.cells[1]])
      
    elif self.mode == 'maxrelbase':
      # Compute release relative to max possible release, subtract base rates.
      factor0 = (np.maximum((rates[self.cells[0]]*self.reldur - self.base_rel[self.cells[0]]), 0.0) / (self.maxrel[self.cells[0]] - self.base_rel[self.cells[0]]))
      factor1 = (np.maximum((rates[self.cells[1]]*self.reldur - self.base_rel[self.cells[1]]), 0.0) / (self.maxrel[self.cells[1]] - self.base_rel[self.cells[1]]))
      
    if self.cell_target == self.cells[0]:
      total_loss = np.mean(np.abs(np.asarray(factor1/factor0))**self.p)
    else:                                 
      total_loss = np.mean(np.abs(np.asarray(factor0/factor1))**self.p)
      
    sample_loss = {'f0': np.mean(np.asarray(factor0)), 'f1': np.mean(np.asarray(factor1)), 'total': total_loss}
      
    # Print?
    if verbose: self.print_loss(sample_loss)
      
    return sample_loss
    
    
############################################################################
############################################################################

class LossOptimizeStimulationMultiParams(LossOptimizeStimulation):
  
  ############################################################################
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
  
  ############################################################################
  def set_base(self, init_rec_data):
    self.base_rates = {cell: np.array([rec_data_i['rate'].mean() for rec_data_i in rec_data['Data']]) for cell, rec_data in init_rec_data.items()}
    self.base_rel = {cell: base_rate*self.reldur for cell, base_rate in self.base_rates.items()}
    
  ############################################################################
  def get_rates(self, rec_data_dict):
    rates = {cell: np.array([rec_data_i.mean() for rec_data_i in rec_data_dict[cell]['rate'][0,:,:].T]) for cell in self.cells}
    return rates
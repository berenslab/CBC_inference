import numpy as np

############################################################################
from delfi.simulator import BaseSimulator
class retsim_Simulator(BaseSimulator):
  
  @copy_ancestor_docstring
  def __init__(self, delfi_optim, **kwargs):
    super().__init__(**kwargs)
    self.delfi_optim = delfi_optim
    
  @copy_ancestor_docstring
  def gen(self, params_list, n_reps=1, pbar=None):
    samples = self.delfi_optim.gen_samples(params_list, n_reps=n_reps)
    if pbar is not None:
      pbar.update(len(params_list))
    return samples
        
############################################################################
from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
class retsim_SummaryStats(BaseSummaryStats):
  
  @copy_ancestor_docstring
  def __init__(self, n_summary=1, seed=None):
    super().__init__(seed=seed)
    self.n_summary = n_summary
    
  @copy_ancestor_docstring
  def calc(self, repetition_list):
    n_reps = len(repetition_list)
    repetition_stats_matrix = np.zeros((n_reps, self.n_summary))
    
    for rep_idx, rep in enumerate(repetition_list):
      repetition_stats_matrix[rep_idx, :] = rep 
    
    return repetition_stats_matrix
        
############################################################################
from delfi.generator.BaseGenerator import BaseGenerator
class retsim_Generator(BaseGenerator):
  
  @copy_ancestor_docstring
  def __init__(self, *args, gen_minibatch, **kwargs):
    super().__init__(*args, **kwargs)
    self.gen_minibatch = gen_minibatch
  
  @copy_ancestor_docstring
  def gen(self, *args, **kwargs):
    if 'minibatch' in kwargs:
      minibatch = kwargs.pop('minibatch')
    else:
      minibatch = self.gen_minibatch
    return super().gen(*args, minibatch=minibatch, **kwargs)
    
  @copy_ancestor_docstring
  def _feedback_proposed_param(self, param):
    return 'accept'

############################################################################    
from delfi import distribution
from TruncatedNormal import TruncatedNormal
def normal2truncated_normal(normal_dist, lower, upper):
    
    if isinstance(normal_dist, distribution.mixture.GaussianMixture.MoG):
        assert normal_dist.ncomp == 1, 'Can not handle MoG with more than 1 comp'
        m, S = normal_dist.calc_mean_and_cov()
    elif isinstance(normal_dist, distribution.Gaussian):
        m = normal_dist.m
        S = normal_dist.S
    else:
        raise NotImplementedError()
        
    return TruncatedNormal(m=m, S=S, lower=lower, upper=upper)
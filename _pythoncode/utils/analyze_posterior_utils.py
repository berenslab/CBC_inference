import numpy as np
import data_utils
import os

##############################################################################################################
def get_samples(posterior, n_samples=200, seed=777, plot=True, prior=None, params=None, plot_lbs=None, plot_ubs=None, plot_opt_x=True):
  posterior.reseed(seed)
  
  o_params = posterior.gen(n_samples)

  if plot:
    assert prior is not None, 'Needed for plotting'
    assert params is not None, 'Needed for plotting'
  
    if plot_lbs is None:
      plot_lbs = np.min(o_params, axis=0)
      plot_ubs = np.max(o_params, axis=0)
    
    import plot_sampling_dists
    
    PP = plot_sampling_dists.SamplingDistPlotter(
      params=params, prior=prior, posterior_list=[posterior],
      lbs=plot_lbs, ubs=plot_ubs
    )
    PP.plot_sampling_dists_1D(
      plot_peak_lines=False, figsize=(12,8), opt_x=plot_opt_x, opt_samples=o_params
    )
    
    return o_params
    
##############################################################################################################
def gen_or_load_samples(optim, opt_params, filename, load):
  if load:
    assert os.path.isfile(filename), 'File does not exist'
    model_output_list = data_utils.load_var(filename)
  else:
    optim.init_rec_data(allow_loading=False, force_loading=True, verbose=True)
    model_output_list = optim.run_parallel(opt_params_list=opt_params, verbose=True)
    data_utils.save_var(model_output_list, filename)

  if load:
    assert len(model_output_list) == opt_params.shape[0], 'Loaded sample size differs from requested'
  
  return model_output_list
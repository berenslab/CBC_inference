from matplotlib import colorbar, colors, cm
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.patches import Ellipse


class SamplingDistPlotter():
  
  ############################################################################
  def __init__(self, params, prior, posterior_list=[], lbs=None, ubs=None, n_x=100):

    self.params = params
    self.prior  = prior
    if isinstance(posterior_list, list):
      self.posterior_list = posterior_list
    else:
      self.posterior_list = [posterior_list]
    
    self.n_x = n_x
    
    if lbs is not None and ubs is not None:
      self.set_bounds(lbs=lbs, ubs=ubs)
 
    self.prior_marginal_1d = None
    self.prior_marginal_2d = None
    self.post_marginal_1d_list = None
    self.post_marginal_2d_list = None

  ############################################################################    
  def set_bounds(self, lbs=None, ubs=None):
    self.bounds = {p_name: (lb, ub) for p_name, lb, ub in zip(self.params.p_names, lbs, ubs)}
  
  ############################################################################
  def compute_values(self, only_1d=False):

    self.opt_x = {}
    self.sim_x = {}
    
    for param, (lb, ub) in self.bounds.items():
      self.opt_x[param] = np.linspace(lb, ub, self.n_x)
      self.sim_x[param] = self.params.opt_param2sim_param(self.opt_x[param], name=param)
    
    # Compute 1D marginals.
    self.prior_marginal_1d = self.compute_all_1d_marginals(self.prior)
    self.post_marginal_1d_list = [self.compute_all_1d_marginals(posterior) for posterior in self.posterior_list]
    
    # Compute maximum of all 1d marginals.  
    self.marginal_1d_max = np.max(
      [np.max(list(self.prior_marginal_1d.values()), axis=None)] +\
      [np.max(list(post_marginal_1d.values()), axis=None) for post_marginal_1d in self.post_marginal_1d_list], axis=None
    )
    print('Maximum 1D value = {:.3g}'.format(self.marginal_1d_max))
    
    if only_1d: return
  
    # Compute 2D marginals for prior.
    self.prior_marginal_2d = self.compute_all_2d_marginals(self.prior)
        
    # Compute 2D marginals for posterior.
    self.post_marginal_2d_list = []
    for posterior in self.posterior_list:
      self.post_marginal_2d_list.append(self.compute_all_2d_marginals(posterior))

    # Compute maximum of all 2d marginals.
    self.prior_marginal_2d_max = self.max_from_2d_marginal(self.prior_marginal_2d)
    self.post_marginal_2d_max  = self.max_from_2d_marginal(self.post_marginal_2d_list)
    self.marginal_2d_max       = np.max([self.prior_marginal_2d_max, self.post_marginal_2d_max])
    print('Maximum 2D prior value = {:.3g}'.format(self.prior_marginal_2d_max))
    print('Maximum 2D post value  = {:.3g}'.format(self.post_marginal_2d_max))

    # Set ticks.
    self.ticks = {param: (np.ceil(self.bounds[param][0]), np.floor(self.bounds[param][1])) for idx, param in enumerate(self.params.p_names)}

  ######################################################################
  def compute_all_1d_marginals(self, dist):
    all_1d_marginals = {param: self.eval_1d_marginal(dist, idx, self.opt_x[param]) for idx, param in enumerate(self.params.p_names)}
    return all_1d_marginals
    
  ######################################################################
  def compute_all_2d_marginals(self, dist):
    all_2d_marginals = {}
    for idx1, param1 in enumerate(self.params.p_names):
      all_2d_marginals[param1] = {}
      for idx2, param2 in enumerate(self.params.p_names):
        if idx2 > idx1:
          all_2d_marginals[param1][param2] =\
            self.eval_2d_marginal(dist, idx1, idx2, self.opt_x[param1], self.opt_x[param2])
    
    return all_2d_marginals
    
  ###########################################################################
  @staticmethod
  def eval_1d_marginal(dist, idx, x):
    return dist.eval(np.atleast_2d(x).T, ii=[idx], log=False)

  ###########################################################################
  @staticmethod
  def eval_2d_marginal(dist, idx1, idx2, x1, x2):
    assert x1.shape == x2.shape
    xx, yy = np.meshgrid(x1, x2)
    assert xx.shape == yy.shape
    pos = np.empty((xx.size,2))
    pos[:, 0] = xx.flatten()
    pos[:, 1] = yy.flatten()
    
    zz = dist.eval(pos, ii=[idx1, idx2], log=False)
    assert xx.size == zz.size, str(xx.size) +' != '+ str(zz.size)
    zz = zz.reshape(xx.shape)
    return (xx, yy, zz)
    
  ######################################################################
  @staticmethod
  def max_from_2d_marginal(marginal_2d_list):
    if not isinstance(marginal_2d_list, list):
      marginal_2d_list = [marginal_2d_list]
  
    maximum = 0
    for marginal_2d in marginal_2d_list:
      for k1, v1 in marginal_2d.items():
          for k2, v2 in v1.items():
              maximum = np.max([np.max(v2), maximum])
    return maximum
  
  ######################################################################
  def cov2cor(self, cov=None):
  
    if cov is None:
      try:
        cov = self.posterior_list[-1].S
      except:
        _, cov = self.posterior_list[-1].calc_mean_and_cov()

    # Compute correlation matrix
    cor = np.zeros(cov.shape)
    for idx1 in range(cor.shape[0]):
      for idx2 in range(cor.shape[1]):
        if idx1 == idx2:
          cor[idx1,idx2] = 1
        else:
          cor[idx1,idx2] = cov[idx1,idx2] / (np.sqrt(cov[idx1,idx1])*np.sqrt(cov[idx2,idx2]))
          
    return cor
  
  ######################################################################
  def plot_correlation(self, cov=None):

    cor = self.cov2cor(cov=cov)
  
    # Plot correlation.
    plt.figure(figsize=(10,10))
    ax = plt.subplot(111)
    plt.imshow(cor, vmin=-1, vmax=1, cmap='coolwarm')
    plt.xticks(np.arange(self.params.p_N))
    plt.yticks(np.arange(self.params.p_N))
    ax.set_xticklabels(self.params.p_names, rotation=90)
    ax.set_yticklabels(self.params.p_names)
    plt.colorbar()
    for idx1, param1 in enumerate(self.params.p_names):
      for idx2, param2 in enumerate(self.params.p_names):
        if idx1 > idx2:
          plt.text(idx1, idx2, "{:.2f}".format(cor[idx1,idx2]),
                  horizontalalignment='center', verticalalignment='center', color='w')
  
  ######################################################################
  def plot_variance(self, cov=None):
  
    if cov is None:
      try:
        cov = self.posterior_list[-1].S
      except:
        _, cov = self.posterior_list[-1].calc_mean_and_cov()
    
    cov_1D = {param: cov[idx, idx] for idx, param in enumerate(self.params.p_names)}
    
    plt.figure(figsize=(10,3))
    ax=plt.subplot(111)
    plt.imshow(np.atleast_2d(list(cov_1D.values())))
    for idx, param in enumerate(self.params.p_names):
      plt.text(idx, 0, "{:.2f}".format(cov_1D[param]),
               horizontalalignment='center', verticalalignment='center', color='w')
    ax.set_xticks(np.arange(0,self.params.p_N))
    ax.set_xticklabels(self.params.p_names, rotation=90)
    ax.set_yticks([])
    plt.colorbar(orientation='horizontal', pad=0.5)
    plt.tight_layout()
    plt.show()
  
  ######################################################################
  def plot_sampling_dists_1D(
    self, opt_x=True, figsize=(12,10), sim_samples=None, opt_samples=None, perc=100,
    plot_peak_lines=True, params=None, post_colors=None
  ):
    
    if self.prior_marginal_1d is None:
      self.compute_values(only_1d=True)
    
    if post_colors is None:
      post_colors = [(i, 0, 0, i) for i in np.linspace(0,1,len(self.posterior_list)+1)[1:]]
    
    if params is None: params = self.params.p_names
    else: params = sorted(params)
    
    if opt_x: x = self.opt_x
    else:     x = self.sim_x
    
    # Plot marginals. 1d and 2d.
    fig = plt.figure(figsize=figsize)
    
    for idx, param in enumerate(params):
    
      ax = plt.subplot(int(np.ceil(np.sqrt(len(params)))), int(np.ceil(np.sqrt(len(params)))), idx+1)
      ax.set_title(param, fontsize=12, verticalalignment='bottom')    
      
      # Set x-limits.
      if opt_x or param not in self.params.p_use_log:
        bounds = (x[param][np.isfinite(x[param])][0], x[param][np.isfinite(x[param])][-1])
      else:
        bounds = (float("{0:.3f}".format(x[param][np.isfinite(x[param])][0])), float("{0:.3f}".format(np.nanpercentile(x[param], perc))))
      
      ax.set_xlim(bounds)
      

      # Plot 1d marginals
      ax.plot(x[param], self.prior_marginal_1d[param], label='prior', c='k', linestyle='--')
      for post_idx, post_marginal_1d in enumerate(self.post_marginal_1d_list):
        label = 'post'
        if len(self.post_marginal_1d_list) > 1: label += str(post_idx+1)
        ax.plot(x[param], post_marginal_1d[param], label=label, c=post_colors[post_idx])

      if plot_peak_lines:
        plt.axvline(x[param][np.argmax(self.prior_marginal_1d[param])], c='k', alpha=0.5)
        for post_idx, post_marginal_1d in enumerate(self.post_marginal_1d_list):
          plt.axvline(x[param][np.argmax(post_marginal_1d[param])], alpha=0.5)
      ax.set_ylim([0, None])
      
      if idx+1==int(np.ceil(np.sqrt(len(params)))): plt.legend(fontsize=8, bbox_to_anchor=(1.2, 1), loc=2, borderaxespad=0)
      

      ax2 = ax.twinx()
      if (sim_samples is not None) or (opt_samples is not None):
        if sim_samples is not None:
          samples_x = sim_samples[param]
          if opt_x: samples_x = self.params.sim_param2opt_param(samples_x, name=param)
        else:
          samples_x = opt_samples[:,idx]
          if not opt_x: samples_x = self.params.opt_param2sim_param(samples_x, name=param)
        
        if samples_x.size < 10:
          for sample_x in samples_x:
            ax2.axvline(sample_x, c='k', linestyle='--', alpha=.5)
        else:
          ax2.hist(samples_x, facecolor='k', alpha=.3, range=bounds, align='mid', bins=21, density=True)
      ax.set_ylim([0, None])
      
    plt.tight_layout()
    
    plt.show()
  
  
  ######################################################################
  def plot_sampling_dists_2D(
    self, opt_x=True, figsize=(12,10), samples=None, perc=100, plot_peak_lines=True,
    only_post_idx=None, params=None, post_colors=None, max_2d=False, max_1d=False
  ):
    
    if self.prior_marginal_2d is None:
      self.compute_values(only_1d=False)
    
    # Maxima for plotting.
    if max_1d is True: max_1d = self.marginal_1d_max
    else:              max_1d = None
    if max_2d is True: max_2d = self.marginal_2d_max
    else:              max_2d = None
    
    if params is None: params = self.params.p_names
    else:              params = [param for param in self.params.p_names if param in params]
    
    if opt_x: x = self.opt_x
    else:     x = self.sim_x
    
    if post_colors is None:
      post_colors = [(i, 0, 0) for i in np.linspace(0.3,1,len(self.posterior_list))]
    
    # Plot marginals. 1d and 2d.
    fig = plt.figure(figsize=figsize)
    
    for col_idx, col_param in enumerate(params):
      for row_idx, row_param in enumerate(params):
    
        ax = plt.subplot2grid((len(params), len(params)), (row_idx, col_idx))
        
        # Get bounds.
        x_bounds = [x[col_param][np.isfinite(x[col_param])][0], x[col_param][np.isfinite(x[col_param])][-1]]
        if not opt_x: x_bounds = [float("{0:.3f}".format(x_bounds[0])), float("{0:.3f}".format(np.nanpercentile(x[col_param], perc)))]
        
        y_bounds = [x[row_param][np.isfinite(x[row_param])][0], x[row_param][np.isfinite(x[row_param])][-1]]
        if not opt_x: y_bounds = y_bounds = [float("{0:.3f}".format(y_bounds[0])), float("{0:.3f}".format(np.nanpercentile(x[row_param], perc)))]
        
        # Set x-ticks.
        ax.set_xlim(x_bounds)
        if row_idx+1 == len(params):
          ax.set_xticks(x_bounds)
          ax.tick_params(axis='x', rotation=90)
        else:
          ax.set_xticks([])
        
        # Plot 1d marginals
        if col_idx == row_idx:
          ax.set_yticks([])     
          
          for post_idx, post_marginal_1d in enumerate(self.post_marginal_1d_list):
            ax.plot(x[col_param], post_marginal_1d[col_param], label = 'post marginals', c=post_colors[post_idx])
          ax.plot(x[col_param], self.prior_marginal_1d[col_param], label='prior', c='k')
          if plot_peak_lines:
            plt.axvline(x[col_param][np.argmax(self.prior_marginal_1d[col_param])], c='k', alpha=0.5)
            for post_idx, post_marginal_1d in enumerate(self.post_marginal_1d_list):
              plt.axvline(x[col_param][np.argmax(post_marginal_1d[col_param])], c=post_colors[post_idx], alpha=0.5)
          ax.set_ylim([0, max_1d])
          
          ax2 = ax.twinx()
          ax2.set_yticks([])
          if samples is not None:
            if samples[col_param].size < 10:
              for sample in samples[col_param]:
                if not opt_x: ax2.axvline(sample, c='k', linestyle='--', alpha=.5)
                else:         ax2.axvline(self.params.sim_param2opt_param(sample, name=col_param), c='k', linestyle='--', alpha=.5)
            else:
              if not opt_x: ax2.hist(samples[col_param], facecolor='k', alpha=.3, range=x_bounds, align='mid', bins=41)
              else:         ax2.hist(self.params.sim_param2opt_param(samples[col_param], name=col_param), facecolor='k', alpha=.3, range=x_bounds, align='mid', bins=41)
          
        else:      
          # Set y-ticks.
          ax.set_ylim(y_bounds)
          if col_idx == 0: ax.set_yticks(y_bounds)
          else:            ax.set_yticks([])
            
          # Plot 2d post marginals    
          if col_idx > row_idx:
            if samples is not None:
              
              if 'params' in samples:
                samples_x = samples['params'][col_param].copy()
                samples_y = samples['params'][row_param].copy()
              elif col_param in samples and row_param in samples:
                samples_x = samples[col_param].copy()
                samples_y = samples[row_param].copy()
                
              if opt_x:
                samples_x = self.params.sim_param2opt_param(samples_x, name=col_param)   
                samples_y = self.params.sim_param2opt_param(samples_y, name=row_param)                   

              plt.scatter(samples_x, samples_y, marker='*', zorder=100, edgecolor='b', facecolor=None, s=1)
            
            if len(self.posterior_list) > 1 and only_post_idx is None:
             
              for post_idx, post_mu in enumerate(self.post_mu_list):
                plt.plot(post_mu[col_idx], post_mu[row_idx], marker='.', c=post_colors[post_idx])
            else:
              if len(self.posterior_list) == 1: post_idx = 0
              else:                             post_idx = only_post_idx
              
              # Get x and y.
              xx = self.post_marginal_2d_list[post_idx][col_param][row_param][0]
              yy = self.post_marginal_2d_list[post_idx][col_param][row_param][1]
              if not opt_x:
                xx = self.params.opt_param2sim_param(xx, col_param)
                yy = self.params.opt_param2sim_param(yy, row_param)
              
              # Plot contour lines of post.
              zz = self.post_marginal_2d_list[post_idx][col_param][row_param][2]
              plt.contour(xx, yy, zz, cmap=cm.gist_heat_r, vmin=0, vmax=max_2d, origin='lower')
              
              if plot_peak_lines:
                plt.axvline(x[col_param][np.argmax(self.post_marginal_1d_list[post_idx][col_param])], c='r', alpha=0.5)
                plt.axhline(x[row_param][np.argmax(self.post_marginal_1d_list[post_idx][row_param])], c='r', alpha=0.5)
          
          # Plot 2d prior marginals
          elif col_idx < row_idx:
            
            # Get x and y.
            xx = self.prior_marginal_2d[col_param][row_param][0]
            yy = self.prior_marginal_2d[col_param][row_param][1]
            if not opt_x:
              xx = self.params.opt_param2sim_param(xx, col_param)
              yy = self.params.opt_param2sim_param(yy, row_param)
            
            # Plot contour lines of prior.
            zz = self.prior_marginal_2d[col_param][row_param][2]
            plt.contour(xx, yy, zz, cmap=cm.Greys, vmin=0, vmax=max_2d, origin='lower')         

            if plot_peak_lines:
              plt.axvline(x[col_param][np.argmax(self.prior_marginal_1d[col_param])], c='k', alpha=0.5)
              plt.axhline(x[row_param][np.argmax(self.prior_marginal_1d[row_param])], c='k', alpha=0.5)

        if row_idx == 0: ax.set_title(col_param, fontsize=12, rotation=90, verticalalignment='bottom')
        if col_idx == 0: ax.set_ylabel(row_param, fontsize=12, rotation=0, horizontalalignment='right')
        
    #plt.tight_layout()
    
    ax = fig.add_axes((1, 0.6, 0.01, 0.3))
    cb = colorbar.ColorbarBase(ax, cmap=cm.gist_heat_r, norm=colors.Normalize(vmin=0, vmax=max_2d), orientation='vertical')
    cb.set_label('p_post')
    
    ax = fig.add_axes((1, 0.15, 0.01, 0.3))
    cb = colorbar.ColorbarBase(ax, cmap=cm.Greys, norm=colors.Normalize(vmin=0, vmax=max_2d), orientation='vertical')
    cb.set_label('p_prior')
    
    plt.show()
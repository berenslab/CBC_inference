from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import interpolation_utils

############################################################################
def set_rcParams():
    sns.set_context('paper')
    sns.set_style('ticks')
    plt.rcParams['axes.linewidth']    = .7
    plt.rcParams['xtick.major.width'] = .5
    plt.rcParams['ytick.major.width'] = .5
    plt.rcParams['xtick.minor.width'] = .5
    plt.rcParams['ytick.minor.width'] = .5
    plt.rcParams['xtick.major.size'] = 2
    plt.rcParams['ytick.major.size'] = 2
    plt.rcParams['xtick.minor.size'] = 1
    plt.rcParams['ytick.minor.size'] = 1
    plt.rcParams['font.size']       = 8
    plt.rcParams['axes.titlesize']  = 10
    plt.rcParams['axes.labelsize']  = 8
    plt.rcParams['legend.fontsize'] = 7
    plt.rcParams['xtick.labelsize'] = 7
    plt.rcParams['ytick.labelsize'] = 7
    plt.rcParams['xtick.major.pad']= 3.5
    plt.rcParams['ytick.major.pad']= 3.5
    plt.rcParams['figure.facecolor']=(0.0,0.0,0.0,0.0)
    plt.rcParams['savefig.facecolor']=(1,1,1,0)
    plt.rcParams["savefig.dpi"] = 1200
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams['axes.unicode_minus']=False
    plt.rcParams['text.usetex'] = False  
    plt.rcParams['figure.dpi'] = 120     # only affects the notebook

############################################################################
def iterate_axes(axs):

  '''Make axes iterable, independent of type.
  
  Paramters:
  
  axs : array or list of matplotlib axes or matplotlib axis
    Axes to apply function to.
  
  '''

  if isinstance(axs, list):
    return axs
  elif isinstance(axs, np.ndarray):
    return axs.flatten()
  else:
    return [axs]

############################################################################
def move_xaxis_outward(axs, scale=5):

  '''Move xaxis outward.
  
  Paramters:
  
  axs : array or list of matplotlib axes.
    Axes to apply function to.
  
  scale : float
    How far xaxis will be moved.
  
  '''

  for ax in iterate_axes(axs):
    ax.spines['bottom'].set_position(('outward', scale))

############################################################################
def adjust_log_tick_padding(axs, pad=2.1):

  ''' Change tick padding for all log scaled axes.
  
  Paramters:
  
  axs : array or list of matplotlib axes.
    Axes to apply function to.
    
  pad : float
    Size of padding.
  
  '''

  for ax in iterate_axes(axs):
    if ax.xaxis.get_scale() == 'log':
      ax.tick_params(axis='x', which='major', pad=pad)
      ax.tick_params(axis='x', which='minor', pad=pad)
        
    if ax.yaxis.get_scale() == 'log':
      ax.tick_params(axis='y', which='major', pad=pad)
      ax.tick_params(axis='y', which='minor', pad=pad)
      
############################################################################
def set_labs(axs, xlabs=None, ylabs=None, titles=None):

  '''Set labels and titles for all given axes.
  
  Parameters:
  
  axs : array or list of matplotlib axes.
    Axes to apply function to.
    
  xlabs, ylabs, titles : str, list of str, or None
    Labels/Titles.
    If single str, will be same for all axes.
    Otherwise should have same length as axes.

  '''
      
  for i, ax in enumerate(iterate_axes(axs)):
    if xlabs is not None:
      if isinstance(xlabs, str): xlab = xlabs
      else:                      xlab = xlabs[i]
      ax.set_xlabel(xlab)
      
    if ylabs is not None:
      if isinstance(ylabs, str): ylab = ylabs
      else:                      ylab = ylabs[i]
      ax.set_ylabel(ylab)
            
    if titles is not None:
      if isinstance(titles, str): title = titles
      else:                       title = titles[i]
      ax.set_title(title)
     
############################################################################
def left2right_ax(ax):
    
    '''Create a twin axis, but remove all duplicate spines.
    
    Paramters:
    
    ax : Matplotlib axis
      Original axis to create twin from.
      
    Returns:
    
    ax : Matplotlib axis
      Twin axis with no duplicate spines.
    '''
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    ax = ax.twinx()
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    return ax
    
############################################################################
def plot_downsampled(ax, x, y, N_plot_points=4000, **kwargs):
    x = np.asarray(x)
    y = np.asarray(y)
    
    if x.size < N_plot_points:
        ax.plot(x, y, **kwargs)
    else:
        new_x = np.linspace(np.min(x), np.max(x), N_plot_points)
        new_y = interpolation_utils.in_ex_polate(x, y, new_x)
        ax.plot(new_x, new_y, **kwargs)
        
############################################################################
def fill_between_downsampled(ax, x, y1, y2, N_plot_points=4000, **kwargs):
    
    x = np.asarray(x)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    
    new_x = np.linspace(x.min(), x.max(), N_plot_points)
    new_y1 = interpolation_utils.in_ex_polate(x, y1, new_x)
    new_y2 = interpolation_utils.in_ex_polate(x, y2, new_x)
    ax.fill_between(new_x, new_y1, new_y2, **kwargs)
    
    
############################################################################
# Source: https://stackoverflow.com/questions/23528477/share-axes-in-matplotlib-for-only-part-of-the-subplots
def set_share_axes(axs, target=None, sharex=False, sharey=False):
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_x_axes.join(target, ax)
        if sharey:
            target._shared_y_axes.join(target, ax)
    # Turn off x tick labels and offset text for all but the bottom row
    if sharex and axs.ndim > 1:
        for ax in axs[:-1,:].flat:
            ax.xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            ax.xaxis.offsetText.set_visible(False)
    # Turn off y tick labels and offset text for all but the left most column
    if sharey and axs.ndim > 1:
        for ax in axs[:,1:].flat:
            ax.yaxis.set_tick_params(which='both', labelleft=False, labelright=False)
            ax.yaxis.offsetText.set_visible(False)
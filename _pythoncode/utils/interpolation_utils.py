import numpy as np
from scipy.interpolate import interp1d, griddata

############################################################################
def in_ex_polate(x_old, y_old, x_new, kind='linear'):
    ''' Interpolate trace, extrapolate with constant if necessary.
    
    Parameters:
    
    x_old : 1d-array
      x-data of trace, e.g. time
    
    y_old : 1d-array
      y-data of trace, e.g. release
      Should have same size as x_old
      
    x_new : 1d-array
      new x-data to inter- or extrapolate y_data on.
      
    kind : str (optional)
      Kind of interpolation, default is 'linear'
      
    Returns:
    
    y_new : 1d-array
      y_data for given x_data
    
    '''


    assert x_old.shape == y_old.shape, str(x_old.shape) + '!=' + str(y_old.shape)

    interp_func = interp1d(x_old, y_old, kind=kind)
    
    interp_idx_bool = np.logical_and((x_new<=np.max(x_old)), (x_new>=np.min(x_old)))
    interp_idx = np.where(interp_idx_bool)[0]
    interp_idx_first = np.min(interp_idx)
    interp_idx_last = np.max(interp_idx)
    
    y_new = np.zeros(x_new.shape)
    if interp_idx.size > 0:
      y_new[interp_idx]= interp_func(x_new[interp_idx])
    y_new[0:interp_idx_first] = y_new[interp_idx_first]
    y_new[interp_idx_last:] = y_new[interp_idx_last]
    
    return y_new
    
############################################################################
def interpolate_xyz2grid(x, y, z, n_pixels=None):
  ''' Interpolate z for given data (x, y and z) onto regular grid.
  Output can be seen as the pixel values of an image.
  
  Parameters:
  
  x, y, z : 1d-array
    Data, should all have the same size.
    
  n_pixels : int or None
    Number of pixels per row and column.
    If None, will use automatic estimate.
    
  Returns:
  
  grid_x, grid_y, grid_z : 2d-arrays
    x, y, and z value for every pixel
  
  '''
  
  # List handling.  
  x = np.asarray(x)
  y = np.asarray(y)
  z = np.asarray(z)
  
  assert x.ndim == 1 and y.ndim == 1 and z.ndim == 1
  assert x.size == y.size and x.size == z.size
  
  if n_pixels is None:
    x_n_pixels = (np.max(x) - np.min(x)) / np.min(np.abs(np.diff(np.unique(x))))
    y_n_pixels = (np.max(y) - np.min(y)) / np.min(np.abs(np.diff(np.unique(y))))
    n_pixels = int(np.ceil(np.max([x_n_pixels, y_n_pixels])))
  
  
  grid_x, grid_y = np.meshgrid(
    np.linspace(np.min(x),np.max(x),n_pixels),
    np.linspace(np.min(y),np.max(y),n_pixels),
  )
  points = [(x_i, y_i) for x_i, y_i in zip(x, y)]
  grid_z = griddata(points, z, (grid_x, grid_y), method='linear')
  
  return grid_x, grid_y, grid_z
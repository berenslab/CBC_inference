import numpy as np

##########################################################################
def continue_time_array(time_array, n_reps=2):

  time_array = np.asarray(time_array)

  idx = time_array.size
  dt = time_array[1] - time_array[0]
  tmax = time_array[-1]
  time_array = np.tile(time_array, n_reps)
  for i_rep in range(n_reps):
    if i_rep > 0:
      time_array[i_rep*idx:] += tmax + dt
    
  return time_array
  
############################################################################
def get_closest_t(t_array, t):
  assert t_array.size > 1
  t_array = np.asarray(t_array)

  if t > np.max(t_array):
    t_idx = t_array.size - 1
  elif t < np.min(t_array):
    t_idx = 0
  else:
    t_idx = np.where(t_array >= t)[0][0]
  
  t = t_array[t_idx]
  
  return t, t_idx
  
############################################################################
def get_and_test_dt(dt=None, t_array=None):

  assert (dt is not None) or (t_array is not None)

  if t_array is not None:
    t_array = np.array(t_array)
    t_array_dts = np.unique(np.abs(np.diff(t_array)))
    
    assert np.isclose(t_array_dts.min(), t_array_dts.max())
    
    if dt is not None:
      assert np.allclose(t_array_dts, dt)
    else:
      dt = np.mean(t_array_dts)
    
  return dt
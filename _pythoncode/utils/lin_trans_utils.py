import numpy as np
from scipy.optimize import minimize

############################################################################
def lin_trans(params, trace):
  ''' Linear transformation of trace.
  First param in params is offset, second is scale.
  '''
  return params[0] + params[1]*trace

############################################################################
def lin_trans_loss(params, trace, target, loss_fun):
  ''' Compute loss after linear transformation.
  '''
  return loss_fun(lin_trans(params, trace), target)
  
############################################################################
def best_lin_trans(trace, target, loss_fun, n_iter=10):
  ''' Compute best linear transformation on trace
  to get close to target, for given loss function.
  
  Parameters:
  
  trace : 1d-array
    Trace to perform linear transformation on.
  
  target : 1d-array, same size as trace
    Target trace for linear transformation
    
  loss_fun : callable
    Loss function to minimize.
    Takes linear transformed trace and target as input.
    
  n_iter : int
    Number of iterations with random intialization.
    
  Returns:
  
  best_trace : 1d-array
    Linerally transformed trace with smallest loss.
  
  best_loss : float
    Loss of best_trace.
  '''

  trace = np.asarray(trace)
  target = np.asarray(target) 
  
  # Summary data.
  trace_rng = np.max(trace) - np.min(trace)
  target_rng = np.max(target) - np.min(target)
  target_mean = np.mean(target)

  # Optimize.
  best_params = np.zeros((n_iter+1,2))
  losses = np.full(n_iter+1, np.inf)
   
  best_params[0,:] = [target_mean, 0]
  losses[0] = lin_trans_loss(
    params=best_params[0,:], trace=trace, target=target, loss_fun=loss_fun
  )
    
  for i in range(n_iter):
    i += 1
    
    if trace_rng == 0:
      trace_factor = 1.0
    else:
      trace_factor = np.random.normal(1,0.1) * target_rng/trace_rng + 1e-9
      
    trace_mean = np.mean(trace_factor*trace)
    trace_offset = np.random.normal(target_mean-trace_mean, np.abs(target_mean-trace_mean))
    
    x_opt = minimize(
      lin_trans_loss, x0=[trace_offset, trace_factor], args=(trace, target, loss_fun),
      bounds=[(None, None), (0, None)],
    )
    
    losses[i] = x_opt.fun
    best_params[i,:] = x_opt.x
  
  if np.any(np.isfinite(losses)):
    best_iter = np.nanargmin(losses)
    best_loss = losses[best_iter]
    best_trace = lin_trans(best_params[best_iter,:], trace)
  else:
    best_trace = trace
    best_loss = loss_fun(trace, target)

  return best_trace, best_loss
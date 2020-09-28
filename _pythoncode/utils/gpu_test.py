import theano
from theano import function, config, shared, tensor
from tensorflow.python.client import device_lib
import numpy as np
import time

def get_available_gpus():
  local_device_protos = device_lib.list_local_devices()
  return [x.name for x in local_device_protos if x.device_type == 'GPU']

def run(verbose=False):
  if verbose: 
    print('Available GPUs:')
    print(get_available_gpus())
    
    print('Thenao version = ' + str(theano.__version__))
    print('Theano device = ' + theano.config.device)
    print('Theano home directory = ' + theano.config.base_compiledir + '. This is where you should create a ".theanorc" file')

  vlen = 10 * 30 * 768  # 10 x #cores x # threads per core
  iters = 500
  
  rng = np.random.RandomState(22)
  x = shared(np.asarray(rng.rand(vlen), config.floatX))
  f = function([], tensor.exp(x))
  if verbose: 
    print(f.maker.fgraph.toposort())
  t0 = time.time()
  for i in range(iters):
      r = f()
  t1 = time.time()
  if verbose:
    print("Looping %d times took %f seconds" % (iters, t1 - t0))
    print("Result is %s" % (r,))
  if np.any([isinstance(x.op, tensor.Elemwise) and ('Gpu' not in type(x.op).__name__) for x in f.maker.fgraph.toposort()]):
    success = False
    if verbose: print('Used the CPU')
  else:
    success = True
    if verbose: print('Used the GPU')
  return success

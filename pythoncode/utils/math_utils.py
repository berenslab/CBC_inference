import numpy as np
    
############################################################################
def normalize(x, lower=0, upper=1):
    assert x.ndim == 1
    assert lower < upper
    
    x = x.copy()
    x -= np.min(x)
    
    if np.max(x) != 0:
        x /= np.max(x)
        
    if lower != 0 or upper != 1:
        x *= upper - lower
        x += lower
    
    return x
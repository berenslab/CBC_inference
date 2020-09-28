import numpy as np

############################################################################ 
from delfi.kernel.BaseKernel import BaseKernel
class BaseNoneKernel(BaseKernel):

    def eval(self, x):
        """Kernel for loss calibration.
        Can deal with Nones by returning NaNs.

        Parameters
        ----------
        x : N x dim
            points at which to evaluate kernel

        Returns
        -------
        weights : N
            normalized to be 1. for x = obs
        """
        assert x.shape[0] >= 1, 'x.shape[0] needs to be >= 1'
        assert x.shape[1] == self.dim, 'x.shape[1] needs to be == self.obs'

        out = np.ones((x.shape[0],))

        for n in range(x.shape[0]):
            if x[n] is None:
                us = np.full(self.obs.shape, np.nan)
            else:
                us = np.dot(self.invH, np.array(x[n] - self.obs).T)

            if self.spherical:
                out[n] = self.normalizer * self.kernel(np.linalg.norm(us))
            else:
                for u in us:
                    out[n] *= self.normalizer * self.kernel(u)

        # check fraction of points accepted
        if self.atleast is not None:
            accepted = out > 0.0
            if sum(accepted) / len(accepted) < self.atleast:
                dists = np.linalg.norm(x - self.obs, axis=1)
                N = int(np.round(x.shape[0] * self.atleast))
                idx = np.argsort(dists)[:N]
                out = np.zeros((x.shape[0],))
                out[idx] = 1.
                return out

        return out

class StrechedGauss(BaseNoneKernel):

    def __init__(self, obs, bandwidth, atleast=None, nan_weight=0.0, max_weight_range=0.0):
        ''' Vector kernel.
        
        Additional Parameters:
            nan_weight : float in [0, 1]
                Weight for NaNs values in kernel evaluation.
        
            max_weight_range : 1d-array withs floats in [0, 1], or float
                Distance to obs, that will result in maximum weight.
                Outside, behaves like a Gaussian.
        '''
        
        self.nan_weight = nan_weight
        
        if isinstance(max_weight_range, (float, int)):
            max_weight_range = np.full(obs.size, max_weight_range, dtype=float)
            
        if isinstance(bandwidth, (float, int)):
            bandwidth = np.full(obs.size, bandwidth, dtype=float)
        
        self.max_weight_range = max_weight_range
    
        super().__init__(obs, bandwidth=bandwidth, spherical=False, atleast=atleast)

        self.vector_kernel = True
    
    def eval(self, x):
        """Kernel for loss calibration.

        Parameters
        ----------
        x : N x dim
            points at which to evaluate kernel

        Returns
        -------
        weights : N
            normalized to be 1. for x = obs
        """
        assert x.shape[0] >= 1, 'x.shape[0] needs to be >= 1'
        assert x.shape[1] == self.dim, 'x.shape[1] needs to be == self.obs'

        out = np.ones((x.shape[0],))

        for n in range(x.shape[0]):
            if x[n] is None:
                us = np.full(self.obs.shape, np.nan)
            else:
                us = np.dot(self.invH, np.array(x[n] - self.obs).T)

            for ii, u in enumerate(us):
                out[n] *= self.normalizer * self.kernel(u, ii)

        # check fraction of points accepted
        if self.atleast is not None:
            accepted = out > 0.0
            if sum(accepted) / len(accepted) < self.atleast:
                dists = np.linalg.norm(x - self.obs, axis=1)
                N = int(np.round(x.shape[0] * self.atleast))
                idx = np.argsort(dists)[:N]
                out = np.zeros((x.shape[0],))
                out[idx] = 1.
                return out

        return out


    def kernel(self, u, ii=0):
        ''' Evaluate kernel at specified value.
        Parameters:
            u : float
                Value to evaluate kernel at.
            ii : int
                Index of max_weight_range
        '''
        if np.isnan(u):
            return self.nan_weight
        elif np.abs(u) <= self.max_weight_range[ii]*self.invH[ii,ii]:
            return 1/np.sqrt(2*np.pi)
        else:
            return 1/np.sqrt(2*np.pi)*np.exp(-0.5*(np.abs(u)-self.max_weight_range[ii]*self.invH[ii,ii])**2)


    def set_bandwidth_ii(self, bandwidth, ii):
        """ Set bandwidth
        
        Parameters
        ----------
        bandwidth : 1d-array
            bandwidth of kernel
            
        ii : int
            Index of bandwidth to update.
        """
        
        self.bandwidth[ii] = bandwidth        
        self.set_bandwidth(self.bandwidth)


    def set_max_weight_range_ii(self, max_weight_range, ii):
        self.max_weight_range[ii] = max_weight_range  
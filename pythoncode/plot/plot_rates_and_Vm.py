import numpy as np
from matplotlib import pyplot as plt

def plot_rates_Vms_iGlus(iGlus, target, rates, Vms, ts_iGlus, ts_rec, final_model_output=None):
    plt.figure(1,(12,10))

    ax1  = plt.subplot(511)
    ax2  = plt.subplot(512)
    ax3  = plt.subplot(513)
    ax4  = plt.subplot(514)
    ax5  = plt.subplot(515)

    axes = [ax1, ax2, ax3, ax4, ax5]
    data = [iGlus, rates, Vms, (rates.T - rates[:,0].T).T, (Vms.T - Vms[:,0].T).T]
    names = ['iGlu', 'rate', 'Vm', 'rate-off', 'Vm-off']
        
    c = 'C0'
    for i, (ax, name) in enumerate(zip(axes, names)):
        
        mu = np.median(data[i], axis=0)
        lb = np.percentile(data[i], axis=0, q=10)
        ub = np.percentile(data[i], axis=0, q=90)   
        
        if i == 0:
            plot_time = ts_iGlus
        else:
            plot_time = ts_rec
        
        ax.set_title(name)
        
        if i == 0:
            ax.plot(plot_time, target, label='target')
        
        ax.plot(plot_time, mu, alpha=0.8, label='median')
        ax.fill_between(plot_time, lb, ub, alpha=0.3, facecolor=c, label='percentiles')
        
        if final_model_output is not None:
            ax.plot(plot_time, final_model_output[name], label='best trace')
        
        ax.legend(loc='lower right')
        
    plt.tight_layout()
import numpy as np
from matplotlib import pyplot as plt

def plot_post_vs_marg_sample_loss(post_loss, post_loss_final_cpl):
    assert post_loss_final_cpl['total'].size == post_loss['total'].size
    
    sort_idx1 = np.argsort(post_loss['total'])
    sort_idx2 = np.argsort(post_loss_final_cpl['total'])

    fig, axs = plt.subplots(2,2,figsize=(12,6),sharex='col')

    for i, sort_idx in enumerate([sort_idx1, sort_idx2]):

        axs[i,0].semilogy(post_loss['total'][sort_idx], label='opt cpl')
        axs[i,0].semilogy(post_loss_final_cpl['total'][sort_idx], label='final cpl')

        idx_g = np.argmax(post_loss_final_cpl['total'][sort_idx] > post_loss['total'][sort_idx])
        axs[i,0].axvline(idx_g)

        axs[i,1].scatter(
            post_loss['total'][sort_idx],
            post_loss_final_cpl['total'][sort_idx],
            s=1, color='red'
        )
        
        axs[i,1].plot(
            [0, post_loss['total'][sort_idx].max()],
            [0, post_loss['total'][sort_idx].max()],
        )

    axs[0,0].legend(loc='upper left')
    plt.show()
    
    
    
    n_smaller_loss = np.sum(post_loss_final_cpl['total'] < post_loss['total'])
    print(f"{(n_smaller_loss/post_loss['total'].size):.2%} are smaller when using final cpl")
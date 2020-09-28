from matplotlib import pyplot as plt
import numpy as np
    
def plot_post_vs_marg(post_model_output_list, marginal_model_output_list):
    loss_names = post_model_output_list[0]['loss'].keys()

    fig, axs = plt.subplots(1, len(loss_names), figsize=(15,4))

    for ax, loss_name in zip(axs, loss_names):
        ax.set_title(loss_name)
        
        post_data_loss = [rec_data['loss'][loss_name] for rec_data in post_model_output_list]
        marg_data_loss = [rec_data['loss'][loss_name] for rec_data in marginal_model_output_list]
        
        ax.hist(post_data_loss, alpha=0.5, label='final', facecolor='C0')
        ax.hist(marg_data_loss, alpha=0.3, label='marginals', facecolor='C1')
        
        ax.legend(loc='upper left', bbox_to_anchor=(0,-0.5))
        ax.axvline(0, c='k', lw=4)
        
        ax.axvline(np.median(np.abs(post_data_loss)), color='C0', ls='--', lw=2)
        ax.axvline(np.median(np.abs(marg_data_loss)), color='C1', ls=':', lw=2)
        
    plt.tight_layout()
    plt.show()
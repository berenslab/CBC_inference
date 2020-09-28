from matplotlib import pyplot as plt
from matplotlib import cm


def plot_eq_rec_data(eq_rec_data_list, parallel_params_list):

    mapper = cm.get_cmap('coolwarm', len(parallel_params_list))

    fig, axs =  plt.subplots(1,4,figsize=(12,3))
    axs[0].set(title='d V', ylabel='mV', xlabel='Time')
    axs[2].set(title='d rate', ylabel='ves/s', xlabel='Time')
    for idx, (rec_data, params) in enumerate(zip(eq_rec_data_list, parallel_params_list)):
        
        axs[0].plot(1e3*(rec_data['BC Vm Soma']-rec_data['BC Vm Soma'].iloc[0]), c=mapper(idx))
        axs[1].plot([idx, idx, idx],
                [rec_data['BC Vm Soma'].min(),
                 rec_data['BC Vm Soma'].iloc[0],
                 rec_data['BC Vm Soma'].max()], '-_', c=mapper(idx))
        
        axs[2].plot(rec_data['rate BC'].mean(axis=1)-rec_data['rate BC'].mean(axis=1).iloc[0], c=mapper(idx))
        axs[3].plot([idx, idx, idx],
                    [rec_data['rate BC'].mean(axis=1).min(),
                     rec_data['rate BC'].mean(axis=1).iloc[0],
                     rec_data['rate BC'].mean(axis=1).max()], '-_', label=params[1]['eqfile'], c=mapper(idx))
    axs[-1].legend(loc='upper left', bbox_to_anchor=(1,1))
    plt.tight_layout()
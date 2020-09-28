import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib import pyplot as plt
import numpy as np

def plot_samples(optim, samples, d_sort_index, cell_target, plot_best_n, plot_worst_n):
    
    from matplotlib import cm
    good_mapper = cm.get_cmap('Greens_r', plot_best_n+3)
    bad_mapper = cm.get_cmap('Reds_r', plot_worst_n+3)
    
    time = optim.rec_data[optim.cells[0].bp_type]['Time'][:samples[optim.cells[0].bp_type]['rate'].shape[1]]

    plt.figure(1,(12,3))
    ax = plt.subplot(121)
    plt.title('Stimulus')
    
    for i in range(plot_best_n):
        stim = optim.stim_generator.create_stimulus(params=samples['params'][d_sort_index[i]])
        plt.plot(optim.stim_time, stim*1e6, label='best '+str(i), c=good_mapper(i))
    
    for i in range(plot_worst_n):
        stim = optim.stim_generator.create_stimulus(params=samples['params'][d_sort_index[-(i+1)]])
        plt.plot(optim.stim_time, stim*1e6, label='worst '+str(i), c=bad_mapper(i))

    plt.ylabel(r'Current ($\mu A$)')

    ax = plt.subplot(122)
    plt.title('rate')
    for cell in optim.cells:
        ls = '-' if cell.bp_type == cell_target else ':'
        for i in range(plot_best_n):
            plt.plot(time, samples[cell.bp_type]['rate'][d_sort_index[i]], ls=ls, color=good_mapper(i))
        for i in range(plot_worst_n): 
            plt.plot(time, samples[cell.bp_type]['rate'][d_sort_index[-(i+1)]], ls=ls, color=bad_mapper(i))

    plt.ylabel(r'rate ($ves/s$)')

    plt.tight_layout()

    plt.show()


def plot_Vm_of_z(rec_time, rec_data_list, ON_cell, OFF_cell, sharey='col', rm_offset=False):
    
    fig, axs = plt.subplots(len(rec_data_list),2, figsize=(12,2*len(rec_data_list)), sharex=True, sharey=sharey)

    cmap = cm.viridis
    norm = Normalize(vmin=0, vmax=80)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    for axs_row, rec_data_stim in zip(axs, rec_data_list):
        for rec_data, ax in zip([rec_data_stim[0], rec_data_stim[5]], axs_row):
            ax.set_title(rec_data['cell'])
            cell = ON_cell if rec_data['cell'] == 'CBC5o' else OFF_cell

            for col in rec_data['Data'].columns:
                if 'Vm' in col:
                    comp = int(col[col.find('(')+1:col.find(')')])
                    color = mapper.to_rgba(cell.comp_data.loc[comp]['comsol_z'])
                    
                    if rm_offset:
                        ax.plot(rec_time, 1e3*(rec_data['Data'].loc[:,col]-rec_data['Data'].loc[0,col]),
                                color=color)
                    else:
                        ax.plot(rec_time, 1e3*rec_data['Data'].loc[:,col],
                                color=color)

            fig.colorbar(mapper, ax=ax, label='z (um)')
    plt.tight_layout()
    plt.show()

    
def plot_synapses(rec_time, rec_data_list, rec_data_list2=None, rec_type='Vm', sharey=True, rm_offset=False):
    
    if rec_data_list2 is None:
        rec_data_list2 = [None] * len(rec_data_list)
    
    fig, axs = plt.subplots(len(rec_data_list), 2, figsize=(12,len(rec_data_list)*1.3),
                            sharex=True, sharey=sharey)
    
    for axs_row, rec_data_stim, rec_data_stim2 in zip(axs, rec_data_list, rec_data_list2):
        
        _plot_synapses_data(
            axs_row, rec_time, rec_data_stim, rec_type=rec_type,
            rm_offset=rm_offset, mean_kw=dict(ls='-', lw=1))
        
        if rec_data_stim2 is not None:
            _plot_synapses_data(
                axs_row, rec_time, rec_data_stim2, rec_type=rec_type,
                rm_offset=rm_offset, mean_kw=dict(ls=':', lw=2))

    plt.tight_layout()
    plt.show()


def _plot_synapses_data(axs_row, rec_time, rec_data_stim, rec_type, rm_offset, mean_kw={}):

    N_param_sets = int(len(rec_data_stim)/2)

    for ax, rec_datas in zip(axs_row, [rec_data_stim[:N_param_sets], rec_data_stim[N_param_sets:]]):
        
        for idx, rec_data in enumerate(rec_datas):
            
            color = 'C'+str(idx)
            
            ax.set_title(rec_data['cell'])

            if rec_type == 'rate':
                cols = [col for col in rec_data['Data'].columns if 'rate' in col]

            elif rec_type == 'Vext':
                cols_V = []
                cols_Vm = []
                for col in rec_data['Data'].columns:
                    if 'rate' in col:
                        node = int(col[col.find('(')+1:col.find(')')])
                        cols_V.append('V (' + str(node) + ')')
                        cols_Vm.append('Vm (' + str(node) + ')')

            else:
                cols = []
                for col in rec_data['Data'].columns:
                    if 'rate' in col:
                        node = int(col[col.find('(')+1:col.find(')')])
                        cols.append(rec_type + ' (' + str(node) + ')')

            if rec_type == 'Ca':
                f = 1e6
            elif 'V' in rec_type:
                f = 1e3
            else:
                f = 1

            if rec_type == 'Vext':
                data = rec_data['Data'].loc[:,cols_V].values - rec_data['Data'].loc[:,cols_Vm].values
            else:
                data = rec_data['Data'].loc[:,cols].values

            mu = np.mean(data*f, axis=1)
            if rm_offset: mu -= mu[0]

            ax.plot(rec_time, mu, **mean_kw, color=color)
            ax.set_ylabel(rec_type)
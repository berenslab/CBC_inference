from matplotlib import pyplot as plt
import time
import numpy as np

def run_and_show_cell(cell, plot=True, verbose=False):
   
    if plot:
        plt.figure(1,(6,6))
        plt.imshow(cell.init_retsim(verbose=False))
        plt.axis('off')
        plt.title(str(cell.n_bc_comps))
        plt.show()

    t0 = time.time()
    rec_data, rec_time, rec_stim = cell.run(update_cell_rec_data=True, verbose=verbose);
    tdur = time.time() - t0
    
    return rec_data, rec_time, rec_stim, tdur
    
def run_and_show_cell(cell, plot=True, verbose=False):
   
    if plot:
        plt.figure(1,(6,6))
        plt.imshow(cell.init_retsim(verbose=False))
        plt.axis('off')
        plt.title(str(cell.n_bc_comps))
        plt.show()

    t0 = time.time()
    rec_data, rec_time, rec_stim = cell.run(update_cell_rec_data=True, verbose=verbose);
    tdur = time.time() - t0
    
    return rec_data, rec_time, rec_stim, tdur
    
def set_and_plot_cpl(cell, cpl_dict):
    cell.update_cpl(**cpl_dict)

    fig, axs = plt.subplots(1,3,figsize=(20,10))
    axs[0].imshow(cell.init_retsim(update=True))
    axs[1].imshow(cell.init_retsim(d=14, update=False))
    axs[2].imshow(cell.init_retsim(d=14, R=True, update=False))
    for ax in axs: ax.axis('off')
    plt.show()
    
def plot_params_vs_tdur(params, tdur_list):
    fig, axs = plt.subplots(1,1,figsize=(12,3))
    
    plt.plot(params, tdur_list, '.-', label='recording')
    plt.plot([params[0], params[-1]], [tdur_list[0], tdur_list[-1]], ':', c='gray', label='linear')
    plt.ylabel('Params')
    plt.ylabel('Time (s)')
    plt.xlim(0,None)
    plt.ylim(0,None)
    plt.legend()
   
    plt.show()
    
def plot_params_vs_time_to_Vm_eq(params, rec_data_list, Vm_eq, plt_Vm_name="Vm Soma"):
    fig, axs = plt.subplots(1,1,figsize=(12,3))
    
    for param, rec_data in zip(params, rec_data_list):
        dist_to_eq = np.abs(rec_data[plt_Vm_name].iloc[0]-Vm_eq)
        plt.semilogy(param, dist_to_eq, 'x')
        print(f'{param} --> error = {dist_to_eq*1e3 : .4f} mV')
        if dist_to_eq == 0.0: plt.axvline(param)

    plt.xlabel('param')
    plt.ylabel('dist to eq [V]')
    plt.show()

def plot_reaching_eq(rec_data, rec_time, rec_stim, plt_Vm_name, ylim=(-0.060,-0.030), Vm_eq=None):

    assert plt_Vm_name in rec_data.columns

    plt.figure(figsize=(12,3))
    plt.subplot(131)
    plt.plot(rec_time, rec_stim)
    plt.subplot(132)
    plt.plot(rec_time, rec_data[plt_Vm_name])
    plt.axhline(rec_data[plt_Vm_name].iloc[-1], c='gray')
    plt.ylim(ylim)
    ax = plt.subplot(133)
    plt.semilogy()
    
    if Vm_eq is None: Vm_eq = rec_data[plt_Vm_name].iloc[-1]
    
    plt.fill_between(rec_time, np.abs(rec_data[plt_Vm_name]-Vm_eq), color='r')
    plt.axhline(0.0001, c='k')
    ax.grid(True)
    plt.tight_layout()
    plt.show()
    
def compare_to_final_cpl(rec_time, rec_data, fc_rec_data):
    
    rate_col_idxs = []
    for idx, col in enumerate(fc_rec_data.columns):
        if 'rate' in col:
            rate_col_idxs.append(idx)
            
    Vm_col_idxs = []
    for idx, col in enumerate(fc_rec_data.columns):
        if 'Vm' in col:
            Vm_col_idxs.append(idx)        
            
    plt.figure(1,(12,8))

    
    MSEs = {}
    
    for idx, (col_idxs, name) in enumerate(zip([rate_col_idxs, Vm_col_idxs], ['rate', 'Vm'])):
        plt.subplot2grid((4,2), (2*idx+0,0))
        plt.title(name)
        for col_idx in col_idxs:
            col = fc_rec_data.columns[col_idx]
            plt.plot(rec_time, fc_rec_data.iloc[:,col_idx], label=col, color='red')
            plt.plot(rec_time, rec_data.iloc[:,col_idx], label=col, color='blue', ls='--')
    
        plt.subplot2grid((4,2), (2*idx+0,1))
        for col_idx in col_idxs:
            col = fc_rec_data.columns[col_idx]
            plt.plot(rec_time, fc_rec_data.iloc[:,col_idx]-rec_data.iloc[:,col_idx], label=col)

        plt.subplot2grid((4,2), (2*idx+1,0), colspan=2)

        mu_rd = np.mean([rec_data.iloc[:,col_idx] for col_idx in col_idxs], axis=0)
        std_rd = np.std([rec_data.iloc[:,col_idx] for col_idx in col_idxs], axis=0)
        
        mu_fc = np.mean([fc_rec_data.iloc[:,col_idx] for col_idx in col_idxs], axis=0)
        std_fc = np.std([fc_rec_data.iloc[:,col_idx] for col_idx in col_idxs], axis=0)

        MSE = np.sum((mu_rd-mu_fc)**2)

        MSEs[name] = MSE
        
        plt.title(f'MSE = {MSE:.4g}')

        plt.plot(rec_time, mu_rd, label='mean rec data', color='red')
        plt.plot(rec_time, mu_fc, ls='--', label='mean target', color='blue')
        
        plt.fill_between(rec_time, mu_rd-std_rd, mu_rd+std_rd, color='red', alpha=0.4, label='_', lw=0)
        plt.fill_between(rec_time, mu_fc-std_fc,mu_fc+std_fc, color='blue', alpha=0.2, label='_', lw=0)
    
    plt.tight_layout()
    plt.show()
    
    return MSEs
    
def run_and_compute_cpl_MSE(cell, test_cpl_dict, fc_rec_data):
    
    set_and_plot_cpl(cell, test_cpl_dict)

    cell.update_cpl(**test_cpl_dict)
    rec_data, rec_time, _, tdur = run_and_show_cell(cell, plot=False)
    print(f'min = {tdur/60}')

    MSEs = compare_to_final_cpl(rec_time, rec_data, fc_rec_data)
    
    return {'cpl_dict': [cpl['c'] for cpl in test_cpl_dict.values()],
            'n_comps': cell.n_bc_comps, 'MSEs': MSEs, 'tdur': tdur}

def plot_cpl2MSE(cpl2MSE):
    fig, axs = plt.subplots(1,3, figsize=(12,3), sharex=True)

    axs[0].set_ylabel('MSE rate')
    axs[1].set_ylabel('MSE Vm')
    axs[2].set_ylabel('n_comps')

    for ax in axs: ax.set_xlabel('run time')

    for cpl2MSE_i in cpl2MSE:
        axs[0].plot(cpl2MSE_i['tdur'], cpl2MSE_i['MSEs']['rate'], marker='x', label=cpl2MSE_i['cpl_dict'])
        axs[1].plot(cpl2MSE_i['tdur'], cpl2MSE_i['MSEs']['Vm'], marker='x', label=cpl2MSE_i['cpl_dict'])
        axs[2].plot(cpl2MSE_i['tdur'], cpl2MSE_i['n_comps'], marker='x', label=cpl2MSE_i['cpl_dict'])
    plt.legend(loc='upper left', bbox_to_anchor=(1,1))
        
    plt.tight_layout()
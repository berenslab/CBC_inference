from scipy.signal import find_peaks as f_peaks
from matplotlib import pyplot as plt
import numpy as np
import math_utils

############################################################################  
def find_peaks(
    trace, prom=None, prom_pos=None, prom_neg=None, c='k', plot=False, verbose=True,
    height=None, height_pos=None, height_neg=None, xlims=[], time=None, label='unknown',
):
    # Get prominence of peaks.
    prom_pos = prom_pos or prom
    prom_neg = prom_neg or prom
    
    height_pos = height_pos or height
    height_neg = height_neg or height
    
    trace = np.asarray(trace)
    time = np.asarray(time)
    
    # Compute.
    peak_indices_pos = f_peaks(trace, prominence=prom_pos, height=height_pos)[0]
    peak_indices_neg = f_peaks(-trace, prominence=prom_neg, height=height_neg)[0]
    peak_indices     = np.append(peak_indices_pos, peak_indices_neg)
    if verbose:
      print('Found ' + str(len(peak_indices_pos)) + ' positive peaks for ' + label)
      print('Found ' + str(len(peak_indices_neg)) + ' negative peaks for ' + label)
    
    # Plot?
    if plot:
      plt.figure(figsize=(12,len(xlims)*1.5))
      for idx, xlim in enumerate(xlims):
        plt.subplot(len(xlims),1,idx+1)
        if idx == 0: plt.title(label)
        plt.xlim(xlim)
            
        plt.plot(time, math_utils.normalize(trace), label='rate', c=c)
        for peak_index in peak_indices_pos:
          plt.axvline(time[peak_index], c=c, alpha=0.5)
        for peak_index in peak_indices_neg:
          plt.axvline(time[peak_index], c=c, linestyle='--', alpha=0.5)
                
      plt.tight_layout()
      plt.show()
        
    return peak_indices, peak_indices_pos, peak_indices_neg
    
############################################################################    
def compare_peaks_in_traces(
    trace_list, time_list, xlims=None, color_list=None, label_list=None, params_dict_list=None, stimulus=None,
    verbose=True, plot=False, plot_single=False, plot_hist=True, dt_max = 0.1, base_trace_i=0,
    figsize=(12,2), outputfile=None, rng=(-0.1, 0.1), mode=1, ignore_rec_times=None,
):
    print('range dt: = ' + str(rng) + ' s')

    # Replace missing parameters.
    label_list = label_list or [str(i) for i in range(len(trace_list))]
    color_list = color_list or ['C' + str(i) for i in range(len(trace_list))]
    
    if not(isinstance(time_list, list)): time_list = [time_list] * len(trace_list)
    xlims = xlims or [(time_list[0].min(), time_list[0].max())]
    
    # Find peaks in traces.
    peaks_idx_dict = {}
    
    for trace_i, trace in enumerate(trace_list):
      try:
        params_dict = params_dict_list[trace_i]
      except:
        params_dict = {}
     
      trace = np.asarray(trace)

      # Find peaks in trace.
      trace_peaks, trace_peaks_pos, trace_peaks_neg = find_peaks(
        trace=trace, time=time_list[trace_i], plot=plot_single,
        verbose=verbose, **params_dict, xlims=xlims, c=color_list[trace_i],
        label=label_list[trace_i],
      )
      # Save to dict.
      peaks_idx_dict[label_list[trace_i]] = {}
      peaks_idx_dict[label_list[trace_i]]['all'] = np.array(trace_peaks)
      peaks_idx_dict[label_list[trace_i]]['pos'] = np.array(trace_peaks_pos)
      peaks_idx_dict[label_list[trace_i]]['neg'] = np.array(trace_peaks_neg)

    # Get time and amp of peaks.
    peaks_t_dict = {}
    peaks_A_dict = {}
    for label_i, label in enumerate(label_list):       
      peaks_t_dict[label] = {}
      peaks_A_dict[label] = {}
      for pn in ['pos', 'neg', 'all']:
        peak_times = np.array(time_list[label_i][peaks_idx_dict[label][pn]])
        peak_amps = np.array(trace_list[label_i][peaks_idx_dict[label][pn]])
        
        if ignore_rec_times is not None:
          valid_idxs = np.ones(peak_times.size, dtype=bool)
          
          for ignore_rec_time in ignore_rec_times:
            invalid_idxs = (peak_times >= ignore_rec_time[0]) & (peak_times <= ignore_rec_time[1])
            valid_idxs = valid_idxs & ~invalid_idxs
          
          peak_times = peak_times[valid_idxs]
          peak_amps = peak_amps[valid_idxs]
          
        peaks_t_dict[label][pn] = peak_times
        peaks_A_dict[label][pn] = peak_amps

    # Get dt of peaks.
    peaks_dt_dict = {}
    peaks_t_base  = peaks_t_dict[label_list[base_trace_i]]
    peaks_A_base  = peaks_A_dict[label_list[base_trace_i]]
    for label_i, label in enumerate(label_list):
      if label_i != base_trace_i:
        
        peaks_dt_dict[label] = {}
    
        # Get time differences.
        if mode == 1:
          for pn in ['pos', 'neg']:
            peaks_dt_dict[label][pn] = []
            for idx0_peak, peak_t_base in enumerate(peaks_t_base[pn]):
              idx1_peak = np.argmin(np.abs(peaks_t_dict[label][pn] - peak_t_base))
              
              # t of peaks.
              t0_peak = peaks_t_base[pn][idx0_peak]
              A0_peak = peaks_A_base[pn][idx0_peak]
              
              t1_peak = peaks_t_dict[label][pn][idx1_peak]
              A1_peak = peaks_A_dict[label][pn][idx1_peak]
              
              # dt of peaks.
              dt_peak = t1_peak - t0_peak
              
              # Check if in range.
              if (dt_peak >= rng[0]) and (dt_peak <= rng[1]):
                peaks_dt_dict[label][pn].append((dt_peak, t0_peak, t1_peak, A0_peak, A1_peak))
        
        # Get time differences.
        elif mode == 2:
          for pn in ['pos', 'neg']:
            peaks_dt_dict[label][pn] = []
            for idx1_peak, peak_t_trace in enumerate(peaks_t_dict[label][pn]):
              idx0_peak = np.argmin(np.abs(peaks_t_base[pn] - peak_t_trace))
              
              # t of peaks.
              t0_peak = peaks_t_base[pn][idx0_peak]
              A0_peak = peaks_A_base[pn][idx0_peak]
              
              t1_peak = peaks_t_dict[label][pn][idx1_peak]
              A1_peak = peaks_A_dict[label][pn][idx1_peak]
              
              # dt of peaks.
              dt_peak = t1_peak - t0_peak
              
              # Check if in range.
              if (dt_peak >= rng[0]) and (dt_peak <= rng[1]):
                peaks_dt_dict[label][pn].append((dt_peak, t0_peak, t1_peak, A0_peak, A1_peak))
        else:
          raise('state mode')
            
        # Concatenate pos and neg.
        peaks_dt_dict[label]['all'] = np.concatenate([peaks_dt_dict[label]['pos'], peaks_dt_dict[label]['neg']])
            
    # Plotting.
            
    # Plot traces in x excerpts.
    if plot:
      plt.figure(figsize=(figsize[0],1+len(xlims)*figsize[1]/2))
      for xlim_i, xlim in enumerate(xlims):
        plt.subplot(len(xlims),1,xlim_i+1)
        plt.xlim(xlim)
        #plt.yticks([])
        
        if stimulus is not None:
          plt.plot(stimulus['Time'], stimulus['Stim'], c='m', alpha=0.3, label='Stim')
        
        for label_i, label in enumerate(label_list):
          time  = time_list[label_i]
          trace = trace_list[label_i]
          color = color_list[label_i]
      
          # Plot trace.
          plt.plot(time, trace, label=label, c=color)
          
          # Plot positive and negative peaks.
          y_heights = {}
          y_heights['neg'] = np.linspace(0.1, -0.2, len(label_list))
          y_heights['pos'] = np.linspace(0.7, 1, len(label_list))
          
          for pn in ['pos', 'neg']:
            if pn == 'pos':
              ls = '-'
              marker = '*'
            else:
              ls = '--'
              marker = 'd'
              
            # Plot peak.
            for peak_t in peaks_t_dict[label][pn]:
                plt.axvline(peak_t, c=color, alpha=0.5, linestyle=ls)
                
            # Plot connection.
            if label in peaks_dt_dict:
              for peak_dt, peak_t_base, peak_t_trace, _, _ in peaks_dt_dict[label][pn]:
                plt.plot([peak_t_base, peak_t_trace], [y_heights[pn][label_i]]*2, ls='-', marker=marker, c=color, alpha=0.5)
        
        if xlim_i == 0:
          plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=11)
      
      if xlim_i +1 == len(xlims): plt.xlabel('Time [s]')
      
      plt.tight_layout()
      if outputfile is not None: plt.savefig(outputfile, dpi=600)
      plt.show()
            
    # Plot histograms?
    if plot_hist:
        
      plt.figure(figsize=figsize)
      
      ax1 = plt.subplot(1,3,1)
      ax2 = plt.subplot(1,3,2)
      ax3 = plt.subplot(1,3,3)
      
      axes = [ax1, ax2, ax3]
      
      base_label = label_list[base_trace_i]
      
      for label_i, label in enumerate(label_list):
        if label_i != base_trace_i:
          # Plot as histograms.
          
          color = color_list[label_i]
          
          for pn_i, pn in enumerate(['pos', 'neg', 'all']):
            peaks = np.array([peak_dt for peak_dt, _, _, _, _ in peaks_dt_dict[label][pn]])
            legend = label + ' - ' + base_label
            ax = axes[pn_i]
            
            # Plot.
            ax.set_title(pn + ' peaks')
            ax.hist(peaks, alpha=0.3, bins=20, range=rng, label=legend, facecolor=color)
            
            peaks_in_range = peaks[(peaks>=rng[0]) & (peaks<=rng[1])]
            
            if peaks_in_range.size > 0:
              ax.axvline(np.median(peaks_in_range), c=color)
            else:
              print('no peaks in range')
          
      ax3.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=11)
      
      for ax in axes:
        ax.set_xlim(rng)
        ax.set_xlabel('dt [s]')
      plt.tight_layout()        
      if outputfile is not None: plt.savefig(outputfile, dpi=600) 
      plt.show()
        
    return peaks_idx_dict, peaks_t_dict, peaks_dt_dict
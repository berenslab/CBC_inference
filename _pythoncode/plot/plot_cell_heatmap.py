import os
import pandas as pd
import numpy as np
from time import sleep
from matplotlib import pyplot as plt
from matplotlib import cm, colors, colorbar, rc
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from PIL import Image, ImageDraw

import time_utils
import data_utils
   
############################################################################
class CellPlotter():
  
  rec_type = None
  colormapping = {}
  morph_data = None
  comp_data = None
  data_draw_morph = None
  data_draw_cone_spheres = None
  
  ############################################################################
  def __init__(
      self, cell, rec_type=None,
      background_color_default = (120, 120, 120, 255),
      cell_color_default= (180, 180, 180, 255),
    ):
    
    self.cell = cell
    
    if (rec_type is None) and (len(self.cell.rec_data.keys()) > 0):
      self.rec_type = list(self.cell.rec_data.keys())[0]
    else:
      self.rec_type = rec_type
      
    assert self.rec_type in self.cell.rec_data.keys(), self.rec_type
  
    if 'comsol_x' not in self.cell.morph_data:
      self.cell.add_comsol_xyz_to_morph_data() 
    self.morph_data = self.cell.morph_data
    
    if 'comsol_x' not in self.cell.comp_data:
      self.cell.add_comsol_xyz_to_comp_data() 
    self.comp_data = self.cell.comp_data
    
    self.get_rec_data_from_cell()
    
    self.set_cell_color_default(cell_color_default)
    self.set_background_color_default(background_color_default)
    
    self.computed_data = False
    
    self.comp_masks = {'im_size': None, 'center': None}
    
    self.node2comp = {}
    for node in self.morph_data['node']:
      self.node2comp[node] = self.cell.node2comp(node)
    
  ############################################################################
  def get_rec_data_from_cell(self):
    ''' Get recording data and time from cell and save.
    Compute Vext, where possible.
    '''
    
    rec_data_raw = self.cell.get_rec_data(rec_type=self.rec_type)
    
    rec_data = {'Vm': {}, 'V': {}, 'Vext': {}, 'rate': {} ,'Ca':{}}
    
    for idx, col in enumerate(rec_data_raw.keys()):
      if 'Vm ' in col:
        comp = int(col[col.find('(')+1:col.find(')')])
        rec_data['Vm'][comp] = rec_data_raw[col].values * 1e3
      elif 'V ' in col:
        comp = int(col[col.find('(')+1:col.find(')')])
        rec_data['V'][comp] = rec_data_raw[col].values * 1e3
      elif 'Ca ' in col:
        comp = int(col[col.find('(')+1:col.find(')')])
        if not np.all(rec_data_raw[col].values == 0):
          rec_data['Ca'][comp] = rec_data_raw[col].values * 1e9
      elif 'rate BC' in col:
        node = int(col[col.find('(')+1:col.find(')')])
        rec_data['rate'][node] = rec_data_raw[col].values
    
    for comp in rec_data['Vm'].keys():
      if comp in rec_data['V'].keys():
        rec_data['Vext'][comp] = rec_data['V'][comp] - rec_data['Vm'][comp]

    self.rec_data = {k: pd.DataFrame(v) for k, v in rec_data.items()}
    
    self.rec_time = self.cell.get_rec_time(rec_type=self.rec_type)
    
  ############################################################################
  def set_cell_color_default(self, cell_color_default):
    ''' Set default cell color.
    '''
    self.cell_color_default = cell_color_default
    
  ############################################################################
  def set_background_color_default(self, background_color_default):
    ''' Set default background color.
    '''
    self.background_color_default = background_color_default
  
  ############################################################################
  def set_colormapping(
      self, plot_type, data=None, symmetric=False,
      y_min=None, y_max=None, cmap=None
    ):
    ''' Prepare colormapping, i.e. how values will translate to color in the plot.
    '''
    
    # Get minimum and maximum Voltage for color range.
    if (y_min is None) and (data is None):
      y_min = np.nanmin(self.rec_data[plot_type].values, axis=None)
    elif (y_min is None) and (data is not None):
      y_min = np.nanmin(data, axis=None)
    assert y_min is not None
    
    if (y_max is None) and (data is None):
      y_max = np.nanmax(self.rec_data[plot_type].values, axis=None)
    elif (y_max is None) and (data is not None):
      y_max = np.nanmax(data, axis=None)
    assert y_max is not None
    
    # Force range > 0
    if y_min == y_max:
      y_min -= 1
      y_max += 1
    y_abs_max = np.max(np.abs([y_min, y_max]))
    
    # Specify the color map.
    if cmap is None:
      if plot_type in ['V', 'Vm']:
        cmap = cm.coolwarm
      elif plot_type == 'Vext':
        cmap = cm.gist_gray
      elif plot_type == 'Ca':
        cmap = cm.Reds
      else:
        cmap = cm.Greens
    
    # Norm and mapper.
    if not(symmetric):
      norm = colors.Normalize(vmin=y_min, vmax=y_max, clip=True)
    else:
      norm = colors.Normalize(vmin=-(y_abs_max), vmax=+(y_abs_max), clip=True)
    
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    # Save.
    self.colormapping[plot_type] = {}
    self.colormapping[plot_type]['cmap'] = cmap
    self.colormapping[plot_type]['norm'] = norm
    self.colormapping[plot_type]['mapper'] = mapper
    self.colormapping[plot_type]['y_min'] = y_min
    self.colormapping[plot_type]['y_max'] = y_max
  
  ############################################################################
  def compute_draw_data(
    self, scale_n=1, inc_cone=False, inc_release=True, release_rad=1,
    x_name='x', y_name='z', flipz=False,
  ):
    ''' Compute geometry of polygons and compartments to plot.
    '''
    
    # Rescale.
    x   = self.morph_data[x_name].values * scale_n
    z   = self.morph_data[y_name].values * scale_n
    dia = self.morph_data['dia'].values  * scale_n
    
    if flipz: z *= -1
    
    # Limits.
    self.x_min = np.min(x)
    self.x_max = np.max(x)
    self.z_min = np.min(z)
    self.z_max = np.max(z)
    
    # Compute geometry of all nodes.
    self.data_draw_morph = {}
    for node_idx in range(self.morph_data.shape[0]):
      node = self.morph_data['node'][node_idx]
      self.data_draw_morph[node] = self.compute_draw_data_node(n_i=node_idx, x=x, z=z, dia=dia)

    self.data_draw_cone_spheres = {}
    if inc_cone:
      for comp, comp_data in self.comp_data.iterrows():
        if comp_data['cell_type'] == 0:
          self.data_draw_cone_spheres[comp] = self.compute_draw_cone_data(comp_data, x_name, y_name, scale_n)

    self.data_draw_release_spheres = {}
    if inc_release:
      for node in self.rec_data['rate'].keys():
        self.data_draw_release_spheres[node] = self.compute_draw_release_data(node, release_rad, x=x, z=z, dia=dia)
 
    self.computed_data = True

  ############################################################################
  def compute_draw_data_node(self, n_i, x, z, dia, plot=False):
    ''' Compute the polygons and sphere for a given node index.
    Note that the node index is not the node.
    '''
    
    node   = self.morph_data['node'][n_i]
    p_list = self.morph_data['p_idx'][n_i]
    c_list = self.morph_data['c_idx'][n_i]
    region = self.morph_data['region'][n_i]    
    
    data_draw_morph_node = {}
    data_draw_morph_node['spheres'] = []
    data_draw_morph_node['polygons'] = []

    for p_i in p_list:
      for c_i in c_list:
      
        # Get points of the nodes.
        ptp = np.array([x[p_i],z[p_i]]).astype(float)
        ptn = np.array([x[n_i],z[n_i]]).astype(float)
        ptc = np.array([x[c_i],z[c_i]]).astype(float)
        
        # Get connecting vectors.
        vec_a = ptn - ptp
        vec_b = ptn - ptc
    
        # Get orthogonal vectors.
        vec_a_o = np.array([-vec_a[1], vec_a[0]])
        vec_b_o = np.array([vec_b[1], -vec_b[0]])
        
        if all(ptn == ptc):
            vec_c = vec_a_o
        else:
            vec_c = (vec_a_o + vec_b_o)
        
        vec_a_norm = np.linalg.norm(vec_a)
        
        vec_a_o_norm = np.linalg.norm(vec_a_o)
        vec_b_o_norm = np.linalg.norm(vec_b_o)
    
        radp = dia[p_i] / 2
        radn = dia[n_i] / 2
        radc = dia[c_i] / 2
        
        assert radp > 0
        assert radn > 0
        assert radc > 0
    
        # Normalize.
        if vec_a_norm != 0:
          vec_a_normalized = vec_a / vec_a_norm
        else:
          vec_a_normalized = vec_a
         
        if vec_a_o_norm != 0:
          vec_a_o /= vec_a_o_norm       
        if vec_b_o_norm != 0:
          vec_b_o /= vec_b_o_norm
        else:
          vec_b_o = vec_a_o
        
        vec_c_norm = np.linalg.norm(vec_c)
        if vec_c_norm != 0:
          vec_c_normalized = vec_c / vec_c_norm
        else:
          vec_c_normalized = vec_a_o
              
        # Calculate coordinates of the compartment.
        if region == 'SOMA':
          x1 = x[n_i] - radn
          x2 = x[n_i] + radn
          
          y1 = z[n_i] - radn
          y2 = z[n_i] + radn
          
          # Append.
          if len(data_draw_morph_node['spheres']) > 0:
            assert np.all(np.array(data_draw_morph_node['spheres']) == np.array([x1, y1, x2, y2]))
          else:
            data_draw_morph_node['spheres'].append(np.array([x1, y1, x2, y2]))
        else:
        
          if self.morph_data['region'][p_i] == 'SOMA':
            xy1 = ptp + vec_a_normalized * np.sqrt(radp**2 - radn**2) + vec_a_o * radn
            xy2 = ptp + vec_a_normalized * np.sqrt(radp**2 - radn**2) - vec_a_o * radn
          
          else:
            xy1 = ptp + 0.5 * vec_a + vec_a_o * 0.5 * (radp + radn)
            xy2 = ptp + 0.5 * vec_a - vec_a_o * 0.5 * (radp + radn)
          
          xy3 = ptn + vec_c_normalized * radn
          xy4 = ptn - vec_c_normalized * radn
          
          xy5 = ptc + 0.5 * vec_b + vec_b_o * 0.5 * (radc + radn)
          xy6 = ptc + 0.5 * vec_b - vec_b_o * 0.5 * (radc + radn)
          
          data_draw_morph_node['polygons'].append([xy1, xy2, xy4, xy6, xy5, xy3])
          
          #data_draw_morph_node['polygons'].append([xy1, xy2, xy3, xy4])
          #data_draw_morph_node['polygons'].append([xy1, xy2, xy4, xy3])
          #data_draw_morph_node['polygons'].append([xy3, xy4, xy5, xy6])
          #data_draw_morph_node['polygons'].append([xy3, xy4, xy6, xy5])
          
    if plot:
      plt.figure(1,(8,8))
      ax = plt.subplot(111)

      ax.axis('equal')

      plt.plot(ptp[0], ptp[1], "rx")
      plt.text(ptp[0], ptp[1], "ptp")

      plt.plot(ptn[0], ptn[1], "rx")
      plt.text(ptn[0], ptn[1], "ptn")

      plt.plot(ptc[0], ptc[1], "rx")
      plt.text(ptc[0], ptc[1], "ptc")

      for i, xy in enumerate([xy1, xy2, xy3, xy4, xy5, xy6]):
          plt.plot(xy[0], xy[1], "k.")
          plt.text(xy[0], xy[1], "xy"+str(i+1))
          
      plt.plot([ptp[0], ptp[0]+vec_a[0]], [ptp[1],ptp[1]+vec_a[1]])
      plt.text(ptp[0]+0.5*vec_a[0], ptp[1]+0.5*vec_a[1], "vec_a")

      plt.plot([ptp[0], ptp[0]+vec_a_o[0]], [ptp[1],ptp[1]+vec_a_o[1]])
      plt.text(ptp[0]+0.5*vec_a_o[0], ptp[1]+0.5*vec_a_o[1], "vec_a_o")

      plt.plot([ptc[0], ptc[0]+vec_b[0]], [ptc[1],ptc[1]+vec_b[1]])
      plt.text(ptc[0]+0.5*vec_b[0], ptc[1]+0.5*vec_b[1], "vec_b")

      plt.plot([ptc[0], ptc[0]+vec_b_o[0]], [ptc[1],ptc[1]+vec_b_o[1]])
      plt.text(ptc[0]+0.5*vec_b_o[0], ptc[1]+0.5*vec_b_o[1], "vec_b_o")

      plt.plot([ptn[0], ptn[0]+vec_c_normalized[0]], [ptn[1],ptn[1]+vec_c_normalized[1]])
      plt.text(ptn[0]+0.5*vec_c_normalized[0], ptn[1]+0.5*vec_c_normalized[1], "vec_c")
          
    return data_draw_morph_node

  ############################################################################
  def compute_draw_cone_data(self, comp_data, x_name, y_name, scale_n):
    ''' Compute cone comp data to draw.
    '''
    x1 = (comp_data[x_name] - comp_data['dia']/2)*scale_n
    x2 = (comp_data[x_name] + comp_data['dia']/2)*scale_n
    
    y1 = (comp_data[y_name] - comp_data['dia']/2)*scale_n
    y2 = (comp_data[y_name] + comp_data['dia']/2)*scale_n
    
    self.update_lims(xs=[x1, x2], zs=[y1, y2])
    
    return np.array([x1, y1, x2, y2])

  ############################################################################
  def compute_draw_release_data(self, node, release_rad, x, z, dia):
    ''' Compute release compartment for node
    '''
    
    node_idxs = np.where(self.morph_data['node'] == node)[0]
    if node_idxs.size == 1:
      node_idx = node_idxs[0]
      
      row = self.morph_data.iloc[node_idx,:]
      
      x1 = x[node_idx] - release_rad
      x2 = x[node_idx] + release_rad
      
      y1 = z[node_idx] - release_rad
      y2 = z[node_idx] + release_rad
      
      self.update_lims(xs=[x1, x2], zs=[y1, y2])

    else:
      raise ValueError
      
    return np.array([x1, y1, x2, y2])
  
  ############################################################################
  def update_lims(self, xs, zs):
    ''' Update plot limits.
    '''
    self.x_min = np.min([self.x_min, np.min(xs)])
    self.x_max = np.max([self.x_max, np.max(xs)])
    self.z_min = np.min([self.z_min, np.min(zs)])
    self.z_max = np.max([self.z_max, np.max(zs)])
  
  ############################################################################
  def get_image_sequence(
      self, t_list, plot_types, extraspace,
      nodes=None, folder=None, filename=None,
      to_array_stack=False
    ):
    if not self.computed_data:
      print('Please compute data first')
      return
    
    for idx, t in enumerate(t_list):
      image = self.draw_heatmap(
        t=t, nodes=nodes, verbose=(idx==0), extraspace=extraspace,
        plot_types=plot_types
      )
      
      if folder is not None:
        assert filename is not None
        assert not to_array_stack, 'Can not use both'
        image.save(folder + '/' + filename + str(idx).zfill(int(np.ceil(np.log10(len(t_list))))) + ".png", "png")
      elif to_array_stack:
        if idx == 0: images = np.empty((len(t_list), ) + np.array(image, dtype=int).shape, dtype=int)
        images[idx] = np.array(image, dtype=int)
      else:
        if idx == 0: images = []
        images.append(image)

      
    return images
  
  ############################################################################
  def get_image_sequences(self, plot_types_list, **kwargs):
    image_sequences = []
    for plot_types in plot_types_list:
      image_sequences.append(self.get_image_sequence(plot_types=plot_types, **kwargs))
    return image_sequences
  
  ############################################################################
  def animate(
      self, data_list, colorbar_list=None, titles=None, cb_labels=None,
      figsize=(12,6), sbnx=2, trace_df=None,
      dt=0.033, filename=None, **kwargs
    ):
    
    dims = data_list[0].shape
    
    from IPython.display import HTML
    import matplotlib.animation as animation
    
    fig, axs, ims, t_line = self.plot_data(
      data_list=data_list, colorbar_list=colorbar_list, titles=titles,
      cb_labels=cb_labels, figsize=figsize, sbnx=sbnx, trace_df=trace_df, **kwargs
    )
    
    def init():
      for im, data in zip(ims, data_list):
        im.set_data(data[0])
         
      return (*ims,)
    
    # animation function. This is called sequentially
    def animate(j):
      for im, data in zip(ims, data_list):
        data_slice = data[j]
        im.set_data(data_slice)
          
      t0 = trace_df.iloc[:,0].min()
      t_rng = trace_df.iloc[:,0].max() - t0
          
      t_line.set_xdata(t0 + t_rng*(j/(data.shape[0]-1)))
          
      return (*ims,)
    
    # call the animator. blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(
      fig, animate, init_func=init, frames=dims[0], interval=dt*1000, blit=True
    )
    
    if filename is not None:
      print('Save to:', filename)
      anim.save(filename)
    
    plt.close()
    return HTML(anim.to_html5_video())

  ############################################################################
  @staticmethod
  def plot_data(
      data_list, colorbar_list=None, titles=None, cb_labels=None,
      figsize=(12,6), sbnx=2, trace_df=None, data_list_idx=0, set_sbny=None,
      cb_width=0.1, abc=None, trace_kw={}, tl_dict={},
    ):
    dims = data_list[0].shape
    
    if colorbar_list is None:
      fig, axs = plt.subplots(1, len(data_list), figsize=figsize, squeeze=False)
      axs = axs.flatten()
      
      assert trace_df is None, 'Not implemented'
    else:
      assert len(colorbar_list) == len(data_list)
      fig = plt.figure(figsize=figsize)
      
      sbny = 5 if trace_df is not None else 1
      
      if set_sbny is not None: sbny = set_sbny
      
      sub_nyx = (sbny,len(data_list)*(sbnx+1))
      
      axs = []
      for i in range(len(data_list)):
        axs.append(plt.subplot2grid(sub_nyx, (0, i*(sbnx+1)), colspan=sbnx, rowspan=sbny-int(trace_df is not None)))
    
      cbaxs = []
      for i in range(len(data_list)):
        cbaxs.append(plt.subplot2grid(sub_nyx, (0, (i+1)*(sbnx+1)-1), rowspan=sbny-int(trace_df is not None)))
    
      if trace_df is not None:
        trace_ax = plt.subplot2grid(sub_nyx, (sbny-1, 0), colspan=sub_nyx[1], rowspan=1)
        trace_ax.axis('off')
        trace_ax.plot(trace_df.iloc[:,0], trace_df.iloc[:,1:], **trace_kw, clip_on=False)
        
        trace_df_min = np.min(trace_df.iloc[:,1:].values, axis=None)
        trace_df_max = np.max(trace_df.iloc[:,1:].values, axis=None)
        trace_df_rng = trace_df_max - trace_df_min
        
        trace_ax.set_ylim(trace_df_min-0.1*trace_df_rng, trace_df_max+0.1*trace_df_rng)
        
        t0 = trace_df.iloc[:,0].min()
        t_rng = trace_df.iloc[:,0].max() - t0          
        t_line = trace_ax.axvline(t0 + t_rng*(data_list_idx/(dims[0]-1)), color='k')
    
    ims = []
    for idx, (ax, data) in enumerate(zip(axs, data_list)):
      ims.append(ax.imshow(data[data_list_idx]))
      
      ax.set_xticks([])
      ax.set_yticks([])
      
      if titles is not None:
        ax.set_title(titles[idx])
      if abc is not None:
        ax.set_title(abc[idx], loc='left', fontweight='bold', ha='right')
    trace_ax.set_title(abc[idx+1], loc='left', fontweight='bold', ha='right')
      
    plt.tight_layout(**tl_dict)
    trace_ax_box = np.array(trace_ax.get_position().bounds)
    trace_ax_box[2] -= axs[0].get_position().bounds[0]
    trace_ax_box[0] = axs[0].get_position().bounds[0]
    trace_ax.set_position(trace_ax_box)
        
        
    if colorbar_list is not None:
      for idx, (ax, cbax) in enumerate(zip(axs, cbaxs)):
    
        colorbar_i = colorbar_list[idx]
        cb = colorbar.ColorbarBase(cbax, cmap=colorbar_i['cmap'], norm=colorbar_i['norm'], orientation='vertical')
        if cb_labels is not None: cb.set_label(cb_labels[idx])
        
        box_im = np.asarray(ax.get_position().bounds)
        box_cb = np.asarray(cbax.get_position().bounds)
        
        box_cb[0] = box_im[0] + box_im[2] + 0.1*box_cb[2]
        box_cb[1] = box_im[1]
        box_cb[2] *= cb_width
        box_cb[3] = box_im[3]
        
        cbax.set_position(box_cb)
        
    return fig, axs, ims, t_line

  ############################################################################
  def draw_heatmap(
      self, plot_types, extraspace, t=None, t_idx=None, nodes=None, verbose=True,
    ):
    ''' Create heatmap image and draw data on it.
    Returns heatmap image.
    '''
  
    if not self.computed_data:
      print('Please compute data first')
      return
  
    assert (t is not None) or (t_idx is not None)
    assert self.morph_data is not None
    
    heatmap = self.create_new_image(extraspace=extraspace)
    center = self.compute_center(extraspace=extraspace)

    # Get time.
    rec_time = self.cell.get_rec_time(rec_type=self.rec_type)
    
    if t_idx is not None:
      assert rec_time[t_idx] == t
    else:
      _, t_idx = time_utils.get_closest_t(t_array=rec_time, t=t)
    
    # Draw cell on image.
    for plot_type in plot_types:
    
      self.compute_masks(
        im_size=heatmap.size, center=center, plot_type=plot_type.split(" ")[0],
        force=False, verbose=False
      )
    
      heatmap = self.draw_cell_on_image(
        image=heatmap, rec_data_t_idx=t_idx,
        center=center, nodes=nodes,
        verbose=verbose, plot_type=plot_type
      )

    return heatmap
    
  ############################################################################
  def compute_center(self, extraspace):
    ''' Compute cell center in image.
    '''
    return (-self.x_min+extraspace, -self.z_min+extraspace)
  
  ############################################################################
  def create_new_image(self, extraspace, color=None):
    ''' Create new image.
    '''
    im_size = np.array([np.max([int((self.x_max-self.x_min)+2*extraspace), 1]),
                        np.max([int((self.z_max-self.z_min)+2*extraspace), 1])])  
    
    if color is None: 
      color = self.background_color_default
    
    return Image.new('RGB', tuple((im_size).astype(int)), color)
    
  ############################################################################
  def draw_cell_on_image(
      self, image, rec_data_t_idx, plot_type, nodes=None, 
      center=(0,0), verbose=True,
    ):
    ''' Draw the cell data on the given image for the given plot_type
    '''
    
    if not self.computed_data:
      print('Please compute data first')
      return

    if plot_type.split(" ")[0] not in self.colormapping:
      print('Please set colormapping first')
      return

    # Remember compartments not found in rec data.
    comps_not_found_in_rec_data = []
    
    # Draw rate compartments?
    if plot_type in ['rate', 'rate with cell']:
      if plot_type == 'rate with cell':
        image, _ = self.draw_compartments_with_color(image=image, nodes=[], rec_data_t_idx=0, center=center, plot_type='Vm')
      self.draw_release_on_image(image, rec_data_t_idx, center)
    
    # Draw node compartments?
    elif plot_type in ['Vext', 'Vm', 'V', 'Ca']:
      image, cnf = self.draw_compartments_with_color(image, nodes, rec_data_t_idx, center, plot_type)
      comps_not_found_in_rec_data += cnf
      
      cnf = self.draw_cone_compartments_with_color(image, rec_data_t_idx, center, plot_type)
      comps_not_found_in_rec_data += cnf
    
    if len(comps_not_found_in_rec_data) > 0 and verbose:
      print(f'The following compartments have no recording data for {plot_type}:')
      print(np.unique(comps_not_found_in_rec_data))
    
    return image    
 
  ############################################################################ 
  def draw_compartments_with_color(self, image, nodes, rec_data_t_idx, center, plot_type):
    ''' Draw cell compartments.   
    '''
    comps_not_found_in_rec_data = []
    drawn_comps = []

    bp_comps = (self.cell.comp_data['cell_type'] == 6) | (self.cell.comp_data['cell_type'] == 10)
    for comp, mask in self.comp_masks[plot_type].items():
      add_node_color = (nodes is None) or (self.cell.comp2node(comp) in nodes)
      if add_node_color:
        if comp in self.rec_data[plot_type].columns:
          color = self.colormapping[plot_type]['mapper'].to_rgba(self.rec_data[plot_type][comp][rec_data_t_idx], bytes=True)
        else:
          comps_not_found_in_rec_data.append(comp)
          color = self.cell_color_default
      else:
        color = self.cell_color_default

      image = Image.composite(Image.new('RGB', image.size, color), image, mask)
          
    return image, comps_not_found_in_rec_data

  ############################################################################
  def draw_release_on_image(self, image, rec_data_t_idx, center):
    ''' Draw the release compartments with imdraw, for given time index and center.
    '''
    
    imdraw = ImageDraw.Draw(image)
    
    for node, data_draw_release_sphere in self.data_draw_release_spheres.items():
      if node in self.rec_data['rate'].columns:
        color = self.colormapping['rate']['mapper'].to_rgba(self.rec_data['rate'][node][rec_data_t_idx], bytes=True)
        imdraw.ellipse([data_draw_release_sphere[i]+center[i%2] for i in [0,1,2,3]], fill=color)
   
  ############################################################################ 
  def draw_cone_compartments_with_color(self, image, rec_data_t_idx, center, plot_type):
    ''' Draw cone spheres with color.
    '''
    comps_not_found_in_rec_data = []
    
    imdraw = ImageDraw.Draw(image)
    
    for comp, data_draw_comp_sphere in self.data_draw_cone_spheres.items():
      if comp in self.rec_data[plot_type].columns:
        color = self.colormapping[plot_type]['mapper'].to_rgba(self.rec_data[plot_type][comp][rec_data_t_idx], bytes=True)
      else:
        comps_not_found_in_rec_data.append(comp)
        color = self.cell_color_default
        
      # Draw sphere.
      imdraw.ellipse([data_draw_comp_sphere[i]+center[i%2] for i in [0,1,2,3]], fill=color)
    
    return comps_not_found_in_rec_data

  ############################################################################
  def compute_masks(self, im_size, center, plot_type, force=False, verbose=False):
    ''' Compute masks to increase plotting speed.
    Values in a mask will always be mapped to the same color.
    e.g. a compartment will have a mask.
    '''
    
    if force:
      recompute_masks = True
    elif (self.comp_masks['im_size'] != im_size) or (self.comp_masks['center'] != center):
      recompute_masks = True
    elif plot_type not in self.comp_masks:
      recompute_masks = True
    else:
      recompute_masks = False
    
    if recompute_masks:
      self.comp_masks['im_size'] = im_size
      self.comp_masks['center'] = center
      self.comp_masks[plot_type] = {}
      
      bp_comps = (self.cell.comp_data['cell_type'] == 6) | (self.cell.comp_data['cell_type'] == 10) 
      
      for comp, comp_data in self.cell.comp_data[bp_comps].iterrows():
        self.comp_masks[plot_type][comp] = self.__compute_comp_mask(comp, plot_type, verbose=verbose)
    
  ############################################################################
  def __compute_comp_mask(self, comp, plot_type, verbose=False):
    ''' Compute the mask of a single compartment.
    '''
    
    center = self.comp_masks['center']
    im_size = self.comp_masks['im_size']
    
    mask = Image.new('1', im_size, 0)
    imdraw = ImageDraw.Draw(mask)
    
    n_spheres = 0
    n_polygons = 0
    n_nodes = 0
    
    for node, data_draw_node in self.data_draw_morph.items():
      if comp == self.node2comp[node]:
        n_nodes += 1
        
        n_spheres += len(data_draw_node['spheres'])
        n_polygons += len(data_draw_node['polygons'])
      
        for sphere in data_draw_node['spheres']:
          imdraw.ellipse([sphere[i]+center[i%2] for i in [0,1,2,3]], fill=1)
        
        for polygon in data_draw_node['polygons']:
          imdraw.polygon([(xy[0]+center[0], xy[1]+center[1]) for xy in polygon], fill=1)
      
    if verbose: print('created mask for comp', comp, ', n_nodes:', n_nodes, ', n_spheres:', n_spheres, ', n_polygons:', n_polygons)
      
    del imdraw
    
    return mask
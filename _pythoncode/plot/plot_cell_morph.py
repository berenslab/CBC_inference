import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

XCONE  = 0
DBP1   = 6
HBP1   = 10
GCA    = 23
GCAOFF = 26
AII    = 13

############################################################################
def plot_2D_cell(
  morph_data=None, comp_data=None, node_list=[], plot_connections=False,
  region_colors={'SOMA':'r','R1':'m','R2':'b','R3':'g','R4':'g','R5':'g','R6':'c','R7':'c','R8':'c'},
  ax=None, inc_cones=False, inc_bc=True, inc_gc=True, plot_node_list_nums=True,
):
  
  assert (morph_data is not None) or (comp_data is not None), 'No data given'
  
  if ax is None:
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
  
  # Plot morphology.
  if morph_data is not None:
    if plot_connections:
      for node_idx in range(morph_data.shape[0]):
        x1 = morph_data['comsol_x'][node_idx]
        y1 = morph_data['comsol_y'][node_idx]
        z1 = morph_data['comsol_z'][node_idx]
        for p_idx in morph_data['p_idx'][node_idx]:
          x2 = morph_data['comsol_x'][p_idx]
          y2 = morph_data['comsol_y'][p_idx]
          z2 = morph_data['comsol_z'][p_idx]
      
          ax1.plot([x1,x2],[y1,y2], c='k', alpha=0.3)
          ax2.plot([x1,x2],[z1,z2], c='k', alpha=0.3)
          ax3.plot([y1,y2],[z1,z2], c='k', alpha=0.3)
            
    xs = morph_data['comsol_x']
    ys = morph_data['comsol_y']
    zs = morph_data['comsol_z']
    cs = [region_colors[morph_data['region'][node_idx]] for node_idx in range(morph_data.shape[0])]
    
    ax1.scatter(xs, ys, c=cs, marker='.', s=morph_data['dia']*100, alpha=0.3)
    ax2.scatter(xs, zs, c=cs, marker='.', s=morph_data['dia']*100, alpha=0.3)
    ax3.scatter(ys, zs, c=cs, marker='.', s=morph_data['dia']*100, alpha=0.3)
  
  # Plot compartments.
  if comp_data is not None:

    assert 'comsol_x' in comp_data, 'Please initialize comsol xyz first to plot comps.'
    
    cone_idxs = (comp_data['cell_type'] == XCONE)
    bp_idxs   = (comp_data['cell_type'] == DBP1) | (comp_data['cell_type'] == HBP1)
    gc_idxs   = (comp_data['cell_type'] == GCA) | (comp_data['cell_type'] == GCAOFF)
    
    for idxs, color in zip([cone_idxs, bp_idxs, gc_idxs], ['red', 'blue', 'orange']):
   
      xs = comp_data['comsol_x'][idxs]
      ys = comp_data['comsol_y'][idxs]
      zs = comp_data['comsol_z'][idxs]
    
      ax1.scatter(xs, ys, c=color, marker='*', s=100, alpha=0.5)
      ax2.scatter(xs, zs, c=color, marker='*', s=100, alpha=0.5)
      ax3.scatter(ys, zs, c=color, marker='*', s=100, alpha=0.5)
    
  for node in node_list:
    node_idxs = np.where(morph_data['node'] == node)[0]
    if node_idxs.size != 1:
      print('{:d} nodes found for node {:d}'.format(node_idxs.size, node))
    else:
      node_idx = node_idxs[0]
      x = morph_data['comsol_x'][node_idx]
      y = morph_data['comsol_y'][node_idx]
      z = morph_data['comsol_z'][node_idx]
      
      ax1.scatter(x, y, marker='.', s=200, c='r', alpha=1)
      if plot_node_list_nums:
        ax1.text(x, y, str(node), zorder=300)
      
      ax2.scatter(x, z, marker='.', s=200, c='r', alpha=1)
      if plot_node_list_nums:
        ax2.text(x, z, str(node), zorder=300)
      
      ax3.scatter(y, z, marker='.', s=200, c='r', alpha=1)
      if plot_node_list_nums:
        ax3.text(y, z, str(node), zorder=300)
      
  ax1.set_xlabel('x')
  ax1.set_ylabel('y')
  
  ax2.set_xlabel('x')
  ax2.set_ylabel('z')
  
  ax3.set_xlabel('y')
  ax3.set_ylabel('z')

  plt.tight_layout()
  plt.show()

############################################################################
def plot_3D_cell(
  morph_data, node_list=[], plot_connections=False,
  region_colors={'SOMA':'r','R1':'m','R2':'b','R3':'g','R4':'g','R5':'g','R6':'c','R7':'c','R8':'c'},
  view=(135, 45), ax=None, plot_node_list_nums=True,
):
  
  if ax is None:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
  
  if plot_connections:
    for node_idx in range(morph_data.shape[0]):
      x1 = morph_data['comsol_x'][node_idx]
      y1 = morph_data['comsol_y'][node_idx]
      z1 = morph_data['comsol_z'][node_idx]
      for p_idx in morph_data['p_idx'][node_idx]:
        x2 = morph_data['comsol_x'][p_idx]
        y2 = morph_data['comsol_y'][p_idx]
        z2 = morph_data['comsol_z'][p_idx]

        ax.plot([x1,x2],[y1,y2],[z1,z2], c='k', alpha=0.3)
          
  xs = morph_data['comsol_x']
  ys = morph_data['comsol_y']
  zs = morph_data['comsol_z']
  cs = [region_colors[morph_data['region'][node_idx]] for node_idx in range(morph_data.shape[0])]
  ax.scatter(xs, ys, zs, c=cs, marker='.', s=morph_data['dia']*100, alpha=0.3)

  for node in node_list:
    node_idxs = np.where(morph_data['node'] == node)[0]
    if node_idxs.size != 1:
      print('{:d} nodes found for node {:d}'.format(node_idxs.size, node))
    else:
      node_idx = node_idxs[0]
      x = morph_data['comsol_x'][node_idx]
      y = morph_data['comsol_y'][node_idx]
      z = morph_data['comsol_z'][node_idx]
      ax.scatter(x, y, z, marker='.', s=200, c='r', alpha=1)
      if plot_node_list_nums:
        ax.text(x, y, z, str(node), zorder=300)
      
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  
  ax.view_init(view[0], view[1])  
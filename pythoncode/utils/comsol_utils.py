import numpy as np
import pandas as pd
import os

############################################################################
def comsol2dataframe(comsol_file):

  assert comsol_file[-4:] == '.csv', 'Please save as csv. '
    
  # Get column names from comments.
  with open(comsol_file, 'r') as f:
    column_names = []
    while ('t=' not in column_names and column_names != ''):
      column_names = f.readline()
  # Clean column names.
  column_names = column_names.split(',')
  column_names[0:3] = ['x', 'y', 'z']
  column_names[-1] = column_names[-1].replace('\n', '')
  # Extract time points.
  time_points = np.zeros(len(column_names)-3)
  for idx, column_name in enumerate(column_names[3:]):
    time_points[idx] = float(column_name[column_name.find('t=')+2:])
  # Reduce column names.
  column_names[3:] = ['t=' + str(time_point) for time_point in time_points]   
  # Load data and transpose.
  data = pd.read_csv(comsol_file, comment='%', header=None)
  data = pd.DataFrame(data.iloc[:,3:].values.T)
  # Update column names.
  data.columns = ['C' + str(i) for i in range(data.shape[1])]
  # Add time.
  data.insert(0, 'Time', time_points)
  
  return data
  
############################################################################
def split_multi_cell_dataframe(Vext_cells_raw, cell, n_expected):
  n_comps = cell.n_bc_comps + cell.n_cone_comps

  time = Vext_cells_raw['Time'].values
  n_cells = int(Vext_cells_raw.shape[1] / n_comps)
  assert n_cells == n_expected, str(n_cells) + ' != ' + str(n_expected)

  assert n_cells*n_comps == Vext_cells_raw.shape[1] - 1, str(n_cells*n_comps) + ' != ' + str(Vext_cells_raw.shape[1] - 1)

  Vext_cells = []
  for cell_i in range(n_cells):
    Vext_cell = Vext_cells_raw.iloc[:,1+cell_i*n_comps:1+(cell_i+1)*n_comps]
    Vext_cell.insert(0, 'Time', time)
    Vext_cells += [Vext_cell]
    
  return Vext_cells

############################################################################
def merge_comsol_comp_files(
    n_comps_per_neuron, prefix, name_list, suffix= '.n', base_folder='Neurons'
  ):
  cells_morph = []
  for name in name_list:
    cell_morph = pd.read_csv(
      os.path.join(base_folder, prefix + name + suffix),
      delim_whitespace=True, names=['x', 'y', 'z']
    )

    assert cell_morph.shape[0] == n_comps_per_neuron
    cells_morph = cells_morph + [cell_morph]

  n_cells = len(cells_morph) 
  cells_morph = pd.concat(cells_morph, axis=0, ignore_index=True)
    
  return cells_morph, n_cells
  
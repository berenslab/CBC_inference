import plot_cell_morph
from shutil import copyfile
import os

comsol2retsim_folder = os.path.join(
  '/gpfs01', 'berens', 'user', 'joesterle', 'berens', 'project_bipolarcells', 'COMSOL2retsim/'
)

############################################################################
def center_xy_region(cell, region, rot=0.0):
  cell.set_comsol_offset(x0=0, y0=0, z0=0)
  cell.set_rot(mxrot=180, myrot=0, mzrot=rot)
  cell.add_comsol_xyz_to_morph_data(verbose=False)

  cell.set_comsol_offset(
    x0=-cell.morph_data["comsol_x"][cell.morph_data["region"] == region].mean(),
    y0=-cell.morph_data["comsol_y"][cell.morph_data["region"] == region].mean(),
  )
  
  cell.add_comsol_xyz_to_comp_data(verbose=False)
    
############################################################################
def create_comsol_comp_file(
    cell, verbose=True, z_min_dist=None, z_soma=None, copy_to_comsol=True,
    rot=0.0, x0=0.0, y0=0.0,
  ):
  assert (z_min_dist is None) or (z_soma is None)
    
  # Correct comps relative to morph.
  cell.auto_correct_comp_data = True
  old_x0 = cell.comsol_x0
  old_y0 = cell.comsol_y0
  cell.set_comsol_offset(x0=0, y0=0, z0=0)
  cell.set_rot(mxrot=180, myrot=0, mzrot=rot)
  
  # Initialize cells.
  cell.init_retsim(verbose=False, print_comps=True, update=True)
  
  # Get z soma.
  if z_min_dist is not None:
    comsol_comps = cell.create_comsol_comp_file(
      comsol_scale_factor=1, verbose=False, inc_cones=True,
      inc_bc=True, inc_gc=False, to_file=False
    )
    z0 = -comsol_comps['z'].min()+z_min_dist
  elif z_soma is not None:
    z0 = z_soma
  else:
    z0 = 0
  
  # Set offset.
  cell.set_comsol_offset(x0=x0+old_x0, y0=y0+old_y0, z0=z0)    
  cell.add_comsol_xyz_to_morph_data(verbose=False)
  cell.add_comsol_xyz_to_comp_data(verbose=False)

  # Create file for COMSOL.
  cell.create_comsol_comp_file(
    comsol_scale_factor=1e-6, verbose=False, inc_cones=True,
    inc_bc=True, inc_gc=False, to_file=True
  )
  
  # Plot compartments.
  if verbose:
    plot_cell_morph.plot_2D_cell(
      morph_data=cell.morph_data, comp_data=cell.comp_data,
      plot_connections=False, inc_cones=True
    )
  
  if copy_to_comsol:
    filename = 'OFF.csv' if cell.is_OFF_bp else 'ON.csv'
    copyfile(
      os.path.join('Neurons', cell.comsol_compfile),
      os.path.join(comsol2retsim_folder, 'comsol_input', 'global', filename)
    );
  else:
    return os.path.join('Neurons', cell.comsol_compfile)
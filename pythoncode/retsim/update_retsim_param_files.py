import os
import data_utils
import numpy as np

def find_and_replace_in_files(
    inputfiles, input_folder, outputfiles, output_folder, params
  ):
  ''' Find params in given files and replace them with new params.
  Save updated files to given output folder.
  '''
  found_params = {}

  data_utils.make_dir(output_folder)
  data_utils.clean_folder(output_folder)

  for inputfile, outputfile in zip(inputfiles, outputfiles):
    
    inputfilepath = os.path.join(input_folder, inputfile)
    assert os.path.isfile(inputfilepath), 'Input file not found: ' + inputfilepath

    with open(inputfilepath, 'r') as f:
      lines = f.readlines()

    # Replace optimized params.
    for p_name, p_value in params.items():
      found_param = False
      
      for idx, line in enumerate(lines):
        if p_name in line:
          found_param = True

          for i in np.flip(np.arange(15)):
            if p_name.ljust(i) in line:
              line = line.replace(p_name.ljust(i), "{:.5g}".format(p_value).ljust(i-1) + " ")

          lines[idx] = line
        
        lines[idx] = lines[idx].replace('\t', ' ')

      if found_param: found_params[p_name] = p_value

    # Write file.
    outputfilepath = os.path.join(output_folder, outputfile)

    print()
    print('saveing updated version of (1) to (2)')
    print('\t(1)', inputfilepath)
    print('\t(2)', outputfilepath)
    print()

    with open(outputfilepath, 'w+') as f:
      for line in lines:
        f.write(line)
        
  print()
  print('The following params were (not) found:')

  max_len_p_name = np.max([len(p_name) for p_name in params.keys()])

  for p_name, p_value in params.items():
    if p_name not in found_params.keys():
      print('\t WARNING: Have not found \t' + p_name.ljust(max_len_p_name) + ' = {:.5g} \t Please replace manually.'.format(p_value))
    else:
      print('\t Found \t' + p_name.ljust(max_len_p_name) + ' = {:.5g}'.format(p_value))
    
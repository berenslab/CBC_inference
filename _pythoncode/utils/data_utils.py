import os
import pickle

############################################################################
def clean_folder(folder, prefix='', exceptions=[], verbose=True, force=False):
  if not(isinstance(exceptions, list)): exceptions = [exceptions]

  files = os.listdir(folder)
  for file in sorted(files):
    if file[0:len(prefix)] == prefix and not(os.path.isdir(os.path.join(folder, file))):
      if (file not in exceptions) and (file not in [exception + '.pkl' for exception in exceptions]):
        file_path = os.path.join(folder, file)
        
        if force:
          os.remove(file_path)
        
        else:
          userinput = input('\t Delete ' + file_path + '? Enter y or yall')
          
          if userinput in ['y', 'yall']:
            os.remove(file_path)
          else:
            print('File', file_path, 'not removed')
            
          if userinput == 'yall':
            force = True

##########################################################################
def is_in_any_file(string, files):
    if not isinstance(files, list): files = [files]
    
    found = False
    for file in files:
      with open(file, 'r') as f:
        lines = f.readlines()
      for line in lines:
        if string in line:
          found = True
          break;
                
    return found

##########################################################################
def make_dir(path):
  if not os.path.exists(path): os.makedirs(path)

############################################################################
def save_var(x, file):
  with open(file, 'wb') as f:
    pickle.dump(x, f, pickle.HIGHEST_PROTOCOL)

############################################################################
def load_var(file):
  with open(file, 'rb') as f:
    x = pickle.load(f)
  return x
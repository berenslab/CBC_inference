import sys
import os

def add2path(path):
  ''' Add path to system path.
  
  path : str
    path to be added to system path.
  
  '''
  if path not in sys.path:
    sys.path = [os.path.abspath(path)] + sys.path
  

def addfolders2path(path, exclude_prefix_list=['.', '_']):

  ''' Add all folders at given path to system path.
  
  path : str
    Path to locate folders.
    path+folder will be added to system path.
  
  exclude_prefix_list : list of chars
    Will not include folders starting with chars in this list.
  '''

  folders = os.listdir(path)
  for folder in folders:
    if folder[0] not in exclude_prefix_list:
      add2path(os.path.join(path, folder))
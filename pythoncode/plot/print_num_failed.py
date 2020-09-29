import numpy as np

def print_num_failed(rec_data_list):
  # Compute number of failed runs.
  failed_list = []
  success_list = []

  total_loss_min = np.inf
  total_loss_max = 0
  total_loss_min_detail = None
  total_loss_max_detail = None

  for i, rec_data in enumerate(rec_data_list):
      if np.any(np.isnan(rec_data['rate'])):
          failed_list.append(i)
      else:
          success_list.append(i)
          
      if rec_data['loss']['total'] > total_loss_max:
          total_loss_max = rec_data['loss']['total']
          total_loss_max_detail = rec_data['loss']
          
      if rec_data['loss']['total'] < total_loss_min:
          total_loss_min = rec_data['loss']['total']
          total_loss_min_detail = rec_data['loss']

  print("N =", len(rec_data_list), "--> {:.0f}% failed".format(100*len(failed_list) / len(rec_data_list)))
  print('\nmax loss:', total_loss_max)
  print_dict(total_loss_max_detail)

  print('\nmin loss:', total_loss_min)
  print_dict(total_loss_min_detail)
  
  return success_list
  
def print_dict(dict):
  for k, v in dict.items():
    print('\t', k, "{:.3g}".format(v))
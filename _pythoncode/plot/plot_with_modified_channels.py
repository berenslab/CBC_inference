from matplotlib import pyplot as plt
import numpy as np
  
############################################################################  
def run_and_plot(optim, final_model_output, multiplier):
    
    sim_params_list, cd_params = get_test_params(final_model_output['params'].copy(), multiplier)

    model_output_list = optim.run_parallel(sim_params_list=sim_params_list, verbose=True)

    plot_comparison_to_final_model_output(
      optim.get_rec_time(), model_output_list, cd_params, final_model_output
    )

############################################################################  
def color_text(txt, color='green'):
    if color == 'green':
        return '\x1b[6;30;42m' + txt + '\x1b[0m'
    elif color == 'red':
        return '\x1b[1;37;41m' + txt + '\x1b[0m'
    else:
        return txt

############################################################################  
def get_test_params(params, multiplier):
  
  cd_params = [param for param in params.keys() if 'cd_' in param]

  # Set every channels to zero and save parameters.
  test_params_list = []

  for param in cd_params:
      test_params_no_cd = params.copy()
      test_params_no_cd[param] *= multiplier
      test_params_list.append(test_params_no_cd)

  for i, param in enumerate(cd_params):
      print(param.ljust(15), end='')
      for j, test_params in enumerate(test_params_list):
          if i == j:
              print(color_text('{:.2f}'.format(test_params[param]).ljust(6)), end=' ')
          else:
              print('{:.2f}'.format(test_params[param]).ljust(6), end=' ')
      print()
      
  return test_params_list, cd_params

############################################################################    
def plot_comparison_to_final_model_output(rec_time, model_output_list, cd_params, final_model_output):

    fig, axs = plt.subplots(len(model_output_list), 2, figsize=(12,1*len(model_output_list)), sharex=True)

    for idx, (model_output, axs_row) in enumerate(zip(model_output_list, axs)):

        if model_output is not None:

            loss_info = "reduced={:.4f}, all={:.4f}".format( model_output['loss']['total'], final_model_output['loss']['total'])
            if model_output['loss']['total'] <= 1.001 * final_model_output['loss']['total']:
                tt_color = 'tomato'
            elif model_output['loss']['total'] <= 1.1*final_model_output['loss']['total']:
                tt_color = 'yellow'
            else:
               tt_color = 'lightgreen'

            axs_row[0].set_title(cd_params[idx][3:] + " - total loss: " + loss_info, backgroundcolor=tt_color)
            axs_row[0].plot(final_model_output['Time'], final_model_output['rate'], label='all channels')
            axs_row[0].plot(rec_time, model_output['rate'], label='_')

            axs_row[1].set_title('Difference')
            
            if np.allclose(final_model_output['Time'], rec_time):
                difference = (model_output['rate'] - final_model_output['rate'])

                if np.all(difference == 0):
                    axs_row[1].text(final_model_output['Time'][0], 0.5, 'No difference', ha='left', va='center')
                    axs_row[1].set_ylim(0,1)
                else:
                    axs_row[1].fill_between(final_model_output['Time'], np.clip(difference, 0, None), color='r', lw=0.0)
                    axs_row[1].fill_between(final_model_output['Time'], np.clip(difference, None, 0), color='b', lw=0.0)
            else:
                axs_row[1].text(final_model_output['Time'][0], 0.5, 'N/A', ha='center', va='center')
                axs_row[1].set_ylim(0,1)

        else:
            print(cd_params[idx])

    axs[0,0].legend()

    plt.tight_layout()
    plt.show()
import numpy as np

def test_if_params_are_used(cell, params, Vm_name, rate_name):
    ''' Change params and check if Vm or rate changes.
    Changes might be small for the given operation range.
    But if they are exactly zero, the parameter is probably not used.
    '''

    all_equal_params = []
    all_close_params = []
    
    for p_name, p_value in params.items():
        
        if isinstance(p_value, (int, float)):
            print('-----------------------------------------------------------')
            print(p_name)

            try:
                if 'off' in p_name:
                    sim_params_list=[{p_name: p_value}, {p_name: p_value-10}, {p_name: p_value+10}]
                elif 'tau' in p_name:
                    sim_params_list=[{p_name: p_value}, {p_name: p_value*3}, {p_name: p_value/3.}]
                elif 'rrp' in p_name:
                    sim_params_list=[{p_name: p_value}, {p_name: 1}, {p_name: 1000}]
                else:
                    sim_params_list=[{p_name: p_value}, {p_name: p_value*5}, {p_name: p_value+1}, {p_name: p_value*0.5}]
            
                test_data_list = cell.run_parallel(
                    plot=True, legend=[p_name],
                    sim_params_list=sim_params_list,
                )

                all_equal = True
                all_close = True

                for test_data_i in test_data_list[1:]:

                    if test_data_i[1].size == test_data_list[0][1].size:
                        if not np.all(test_data_list[0][0][Vm_name].values ==\
                                      test_data_i[0][Vm_name].values):
                            all_equal = False
                        elif not np.all(test_data_list[0][0][rate_name].values ==\
                                        test_data_i[0][rate_name].values):
                            all_equal = False

                        if not np.allclose(test_data_list[0][0][Vm_name].values,
                                           test_data_i[0][Vm_name].values):
                            all_close = False
                        elif not np.allclose(test_data_list[0][0][rate_name].values,
                                             test_data_i[0][rate_name].values):
                            all_close = False
                    else:
                        all_equal = False
                        all_close = False

                if all_equal:
                    all_equal_params.append(p_name)
                    print('warning: all exactly equal')

                elif all_close:
                    all_close_params.append(p_name)
                    print('warning: all close')
                    
            except KeyboardInterrupt:
                return
            except:
                print('Error')
                
    return all_equal_params, all_close_params
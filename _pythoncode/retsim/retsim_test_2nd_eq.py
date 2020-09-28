class test_class():
    def __init__(self):
        pass
    
    @staticmethod
    def is_success(rec_data):
        
        if rec_data is not None:
            if rec_data['BC Vm Soma'].values[0] > -0.030:
                return False
            elif rec_data['BC Vm Soma'].values[-1] > -0.030:
                return False
            else:
                return True
        else:
            return False
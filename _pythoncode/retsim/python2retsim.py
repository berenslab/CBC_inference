import subprocess
import pandas as pd
import numpy as np
from io import StringIO, BytesIO
import re
from PIL import Image

global __retsim_dir
global __vid_dir

def set_dir(retsim_directory, vid_directory=None):
    global __retsim_dir
    global __vid_dir
    __retsim_dir=retsim_directory
    if not vid_directory is None:
        __vid_dir=vid_directory

def retsim(
        retsim_dir=None, vid_dir=None, expt=None, print_labels=False, filename=None,
        errfilename=None, d=None, R=False, vid_args=None, pov_args=None, pov_fn='temp', 
        pov_in_fn=None, pov_image='temp.png', im_size=1000, print_err=0, timeout=99999,
        **kwargs
    ):
    
    if retsim_dir is None:
        if '__retsim_dir' in globals():
            retsim_dir=__retsim_dir
        else:
            print('Retsim directory has to be specified')
            return
    
    if vid_dir is None:
        if '__vid_dir' in globals():
            vid_dir=__vid_dir
        
    arg_list=['./retsim']
    if not expt: #print help if no experiment name is given
        p=subprocess.Popen(arg_list,cwd=retsim_dir,stderr=subprocess.PIPE,universal_newlines=True)
        output=p.communicate(timeout=timeout)
        print(output[1])
        return
        #return output[1] #maybe replace this by own help text
    
    arg_list+=['--expt',expt]
    for key in kwargs.keys():
        if key != 'other_params':
            arg_list+=['--'+key,str(kwargs[key])]
        else:
            for other_key in kwargs['other_params']:
                other_v = kwargs['other_params'][other_key]
                if isinstance(other_v, np.ndarray):
                    if other_v.size == 1:
                        other_v = other_v[0]
                    else:
                        print('Other params must be scalars of arrays with only one value.')
                        return
                   
                arg_list+=['--'+other_key,str(other_v)]
    
    if R: #render povray image
        if not d:
            d=1
        arg_list+=['-d',str(d),'-R','-1',pov_fn+'.pov']
        subprocess.call(arg_list,cwd=retsim_dir)
        pov_arg_list=['povray']
        if not pov_in_fn:
            pov_in_fn=pov_fn
        if pov_args:
            pov_arg_list+=pov_args+['+i'+pov_in_fn+'.pov','+o'+pov_image]
        else:
            pov_arg_list+=['+h'+str(im_size),'+w'+str(im_size),'+i'+pov_in_fn+'.pov','+o'+pov_image]
        print(pov_arg_list)
        subprocess.call(pov_arg_list,cwd=retsim_dir)
        im=Image.open(retsim_dir+pov_image)
        return im, None
    
    elif d: #return image of the model
        arg_list+=['-d',str(d),'-v']
        if not vid_dir:
            vid_dir=retsim_dir.split('models')[0]+'bin/'
        p=subprocess.Popen(arg_list,cwd=retsim_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=False)
        output=p.communicate(timeout=timeout)
        err_out=output[1].decode('utf8')
        vid_arg_list=[vid_dir+'vid','-w',str(im_size)]
        if vid_args:
            vid_arg_list+=vid_args
        q=subprocess.Popen(vid_arg_list+['-c'],cwd=retsim_dir,stdin=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=False)
        output2=q.communicate(input=output[0])
        im=Image.open(BytesIO(output2[0]))

        if print_err:
            return [im,err_out]
        return im
    
    elif filename: #write output to file
        arg_list+=['-1',filename]
        if errfilename:
            arg_list+=['-2',errfilename]
        subprocess.call(arg_list,cwd=retsim_dir)
        return
    
    else: #standard case: return data as pandas dataframe
        p=subprocess.Popen(arg_list,cwd=retsim_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True)
        try:
          output=p.communicate(timeout=timeout)
          
          err_out = output[1]
          data=output[0]
          data=re.sub(' +',' ',data)
          data=re.sub(' \n','\n',data)
          
          comment=''
          for line in data.splitlines():
            if line[0:2]=='#n':
              comment+=line+'\n'
          data=pd.read_csv(StringIO(data),comment='#',delimiter=' ',header=None)

        except:
          try:
            err_out = output[1]
          except:
            err_out = 'Error'
          data=[]
          comment = ''

        labels=comment.replace('#n ','').replace('"','').split('\n')[:-2]
        if (print_labels==1)&(print_err==1):
          return [data,labels,err_out]
        elif print_labels==1:
          return [data,labels]
        elif print_err==1:
            return [data,err_out]
        return data


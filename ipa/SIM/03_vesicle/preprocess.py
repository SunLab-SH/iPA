'''
Extract info from Imaris data files

'''


#%%
import numpy as np
import pandas as pd
import glob, time
import re


#%%
# file_dir = 'G:\\BING\\Fluorescence\\results\\'
file_dir = 'I:\\Bing\\fluorescence\\Results\\'

conditions = ['2.8-30', '16.7-5', '16.7-30']

#%%
# 3D coordinates of ISG, um

for i in conditions:
    tag = file_dir + i + '\\'
    file_dir_list = glob.glob(tag + '*')
    for file_dir_dir in file_dir_list:
        
        file_name = glob.glob(file_dir_dir + '\\' + '*_Position.csv')
        # print(file_name[0])
        position_file = pd.read_csv(file_name[0], header=2) 
        xyz = np.array((position_file.iloc[:,0], position_file.iloc[:,1], position_file.iloc[:,2], position_file.iloc[:,7]-1e9), dtype= "float16")
        
        ret = re.match(r'.*[P,t](.*)_Statistics\\(.*).csv', str(file_name[0]))
        label = ret.group(2)
        print(ret.group(2))

        np.save(file=f'{file_dir}\\coordinates\\{label}.npy', 
                arr=xyz)
        
# %%     
# 3D velocity of ISG, um/s
   
for i in conditions:
    tag = file_dir + i + '\\'
    file_dir_list = glob.glob(tag + '*')
    for file_dir_dir in file_dir_list:

        file_name = glob.glob(file_dir_dir + '\\' + '*_Velocity.csv')
        # print(file_name[0])
        position_file = pd.read_csv(file_name[0], header=2) 
        vel = np.array((position_file.iloc[:,0], position_file.iloc[:,1], 
                        position_file.iloc[:,2], position_file.iloc[:,7]-1e9), dtype= "float16")  # dt = 6 s
        
        ret = re.match(r'.*[P,t](.*)_Statistics\\(.*).csv', str(file_name[0]))
        label = ret.group(2)
        print(label)

        np.save(file=f'{file_dir}\\velocity\\{label}_3dvelocity.npy', 
                arr=vel)


# %%

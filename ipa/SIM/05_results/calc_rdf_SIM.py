#!/usr/bin/python

#%%
# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
from operator import sub
from numpy import zeros, ones, arange, asarray, concatenate
import time
import cv2
#%%
file_dir = 'I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\04_rdf\\' 
file_dir2 = 'I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\00_index\\' 
N_shell = 10

file_list=[]
with open(f'../file_list_MT.txt','r') as f:                           # Modified
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
    print(file_name)

    start_time = time.time() # Measure time elapsed in Python

    '''
    Read data file
    '''

    
    '''
    calculate the rdf profile
    '''

#########################################
    k1 = np.loadtxt(f'{file_dir2}{file_name}_bin2_ne_index.xvg')
    xs = [item[0] for item in k1]
    ys = [item[1] for item in k1]
    zs = [item[2] for item in k1]
    k2 = np.loadtxt(f'{file_dir2}{file_name}_bin2_pm_index.xvg')
    xs2 = [item[0] for item in k2]
    ys2 = [item[1] for item in k2]
    zs2 = [item[2] for item in k2]

    pm_edge = np.zeros((64, 800, 800))   # Modified
    ne_edge = np.zeros((64, 800, 800)) 

    for i in range(0, len(k2)):
        pm_edge[int(xs2[i]),int(ys2[i]),int(zs2[i])] = 1
    for i in range(0, len(k1)):
        ne_edge[int(xs[i]),int(ys[i]),int(zs[i])] = 1

    cell_total_vol = 0
    for z in range(64):
        # Create an empty mask with the same shape as the contour line array
        mask = np.zeros_like(pm_edge[z], dtype=np.uint8)
        # Find contours from the binary contour line array
        contours, _ = cv2.findContours(pm_edge[z].astype('uint8'), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        # Fill the interior of each contour with white color in the mask
        cv2.drawContours(mask, contours, -1, (255), cv2.FILLED, offset=(0, 0))
        mask_array = np.asarray(mask/255)

        mask_ne = np.zeros_like(ne_edge[z], dtype=np.uint8)
        contours_ne, _ = cv2.findContours(ne_edge[z].astype('uint8'), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(mask_ne, contours_ne, -1, (255), cv2.FILLED, offset=(0, 0))
        mask_ne_array = np.asarray(mask_ne/255)

        mask_vol = np.sum(mask_array) - np.sum(mask_ne_array)
        cell_total_vol += mask_vol
#########################################

    # with open(f'{file_dir}{file_name}_bin2_rdf.xvg') as fp:
    with open(f'{file_dir}{file_name}_bin2_rdf_mt.xvg') as fp:
        k1 = fp.readlines()

        
    N_isg = len(k1)/N_shell
    print(len(k1), N_isg)


    isg_in_shell = []
    for i in range(N_shell):
        isg = N_isg - k1[int(N_isg*1*i) : int(N_isg*1*(i+1))].count('\n')
        isg_in_shell.append(isg)
    
    # isgout=open(f'{file_dir}{file_name}_bin2_isgrdf_ne-pm.xvg', 'w') # output
    isgout=open(f'{file_dir}{file_name}_bin2_mtrdf_ne-pm.xvg', 'w') # output

    isg = np.zeros(N_shell-1)
    #print(k1[1,0])

    print(f'isg total number: {isg_in_shell[-1]}')

    for i in range(2, N_shell-1):
        isg[i-2] = abs(isg_in_shell[i]-isg_in_shell[i-1])/isg_in_shell[-1]#/cell_total_vol

    isg[N_shell-3] = abs(isg_in_shell[0]-isg_in_shell[N_shell-2])/isg_in_shell[-1]#/cell_total_vol
    isg[N_shell-2] = abs(isg_in_shell[N_shell-1]-isg_in_shell[0])/isg_in_shell[-1]#/cell_total_vol
    print(f'mean: {np.mean(isg)}')

    np.savetxt(isgout, isg, fmt='%.4f')
    isgout.close()

print("--- %s seconds ---" % (time.time() - start_time))


# %%

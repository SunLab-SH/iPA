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
import os, sys, glob
import re

# os.chdir("I:\\Bing\\fluorescence\\Results\\SIM")   

# os.getcwd() 

# define functions
'''generate n equally distributed points along a 3D vector from point p1 to point p2.'''

def pointsalongvector(ne, pm, n, center, edge = 100):
    points=[]
    v1 = np.linalg.norm(list(map(sub, ne, center)))
    v2 = np.linalg.norm(list(map(sub, pm, center)))
    if v2-edge/31.3 > v1: # pm > ne
        # p2 = [center[0] + v2[0]*((np.linalg.norm(v2)-edge/31.3) / np.linalg.norm(v1)),
        #       center[1] + v2[1]*((np.linalg.norm(v2)-edge/31.3) / np.linalg.norm(v1)),
        #       center[2] + v2[2]*((np.linalg.norm(v2)-edge/31.3) / np.linalg.norm(v1))]
        p1 = ne
        p2 = [ne[0]+(pm[0]-ne[0])*(v2-v1-edge/31.3)/(v2-v1), \
              ne[1]+(pm[1]-ne[1])*(v2-v1-edge/31.3)/(v2-v1), \
              ne[2]+(pm[2]-ne[2])*(v2-v1-edge/31.3)/(v2-v1)]
    elif v2 > v1:
        p1 = ne
        p2 = pm
    else:  # pm <= ne, ne should be cut
        p1 = [pm[0]-1, pm[1]-1, pm[2]]
        p2 = pm
    points.append(p2) # save the pm edge before cut 100 nm

    for i in range(0,n+1):
        points.append([p1[0]+(p2[0]-p1[0])*i/n,p1[1]+(p2[1]-p1[1])*i/n,p1[2]+(p2[2]-p1[2])*i/n])
    return points

def cut_docked(ne,pmpoint,center, edge = 100):
     for i in range(0,len(ne)):
        v1 = list(map(sub, ne[i], center))
        v2 = list(map(sub, pmpoint, center))
        if np.linalg.norm(v2)-edge/31.3 > np.linalg.norm(v1): # pixel
             pmpoint_new = center + v2*((np.linalg.norm(v2)-edge/31.3) / np.linalg.norm(v1))


# seqs = ["2.8-30_P3S3-1_PM_NE_mask", "2.8-30_P3S3-2_PM_NE_mask", 
#         "2.8-30_P17S1-1_PM_NE_mask", "16.7-5_P4_PM_NE_mask", 
#         "16.7-30_P14S1-2_PM_NE_mask", "16.7-30_P22_PM_NE_mask"]

start_time = time.time() # Measure time elapsed in Python
nslice = 8 # 16
num_file = 7
#%%

file_dir = 'I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\00_index\\'  # Modified
file_dir2 = 'I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\01_findpair\\'  # Modified
output_dir = 'I:\\Bing\\fluorescence\\3D\\20230320MT-ISG\\02_slice\\'  # Modified

file_list=[]
with open(f'../file_list_MT.txt','r') as f:                           # Modified
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list:
    print(file_name)

    # Read data file
    p = np.loadtxt(f'{file_dir}{file_name}_bin2_ne_index.xvg')
    center = [float(p[0][0]),float(p[0][1]),float(p[0][2])] # Get center, ne, and pm

    k = {}
    for i in range(num_file):
        # print(i)
        k[i] = np.loadtxt(f'{file_dir2}{file_name}_bin2_pair_ne-pm_shell{i}.xvg')
        # k[1] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_2.xvg')
        # k[2] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_3.xvg')
        # k[3] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_4.xvg')
        # k[4] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_5.xvg')
        # k[5] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_6.xvg')
        # k[6] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_7.xvg')
        # k[7] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_8.xvg')

    # 16 slices
    points=[]

    for i in range(0,num_file):
        for j in range(0,len(k[i])):
            ne = [k[i][j,0],k[i][j,1],k[i][j,2]]
            pm = [k[i][j,3],k[i][j,4],k[i][j,5]]
            points.append(pointsalongvector(ne, pm, nslice, center)[:])

    shapepoints = np.reshape(points, (len(points),nslice+1+1,3)) # Modified
    for i in range(0,nslice+1+1): # it is the NE layer when i = 1  
        #print(i)
        toclean = list(shapepoints[:,i,:]) # to remove repeating sublists
        cleaned = [list(i) for i in set(map(tuple, toclean))]
        # calculate volume of iregular made of n points in a 3D space
        with open(f'{output_dir}{file_name}_bin2_slice_shell{i}.xvg', 'w') as outfile:
            np.savetxt(outfile, cleaned, fmt='%.4f')

# delete repeating points
print("--- %s seconds ---" % (time.time() - start_time))

# %%

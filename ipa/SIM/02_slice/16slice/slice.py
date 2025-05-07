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

os.chdir("I:\\Bing\\fluorescence\\Results")   

os.getcwd() 

# define functions
'''generate n equally distributed points along a 3D vector from point p1 to point p2.'''
def pointsalongvector(p1, p2, n):
    points=[]
    for i in range(0,n+1):
        points.append([p1[0]+(p2[0]-p1[0])*i/n,p1[1]+(p2[1]-p1[1])*i/n,p1[2]+(p2[2]-p1[2])*i/n])
    return points

# seqs = ["2.8-30_P3S3-1_PM_NE_mask", "2.8-30_P3S3-2_PM_NE_mask", 
#         "2.8-30_P17S1-1_PM_NE_mask", "16.7-5_P4_PM_NE_mask", 
#         "16.7-30_P14S1-2_PM_NE_mask", "16.7-30_P22_PM_NE_mask"]
start_time = time.time() # Measure time elapsed in Python
nslice=16
#%%
# Read data file
k = {}

# root_dir = 'F:\\Bing\\Features\\Results\\02_findpair\\outputs\\'
timeslice = 8

seq_list = glob.glob('I:\\Bing\\fluorescence\\3+1D\\*tif')
for seq_file in seq_list:
    seq = seq_file[26:-4]

    print(seq)

    for i in range(nslice):
        k[i] = np.loadtxt(f'02_findpair\\outputs_t8\\{seq}_PM_NE_mask_pair_ne-pm_shell{i}_t{timeslice}.xvg')
        # k[1] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_2.xvg')
        # k[2] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_3.xvg')
        # k[3] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_4.xvg')
        # k[4] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_5.xvg')
        # k[5] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_6.xvg')
        # k[6] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_7.xvg')
        # k[7] = np.loadtxt('../../findpair/output/'+ str(seq) +'_pair_ne-pm_8.xvg')

    # 16 slices
    points=[]

    for i in range(0,nslice):
        for j in range(0,len(k[i])):
            ne = [k[i][j,0],k[i][j,1],k[i][j,2]]
            pm = [k[i][j,3],k[i][j,4],k[i][j,5]]
            points.append(pointsalongvector(ne, pm, nslice)[:])

    shapepoints = np.reshape(points, (len(points),nslice+1,3))
    for i in range(0,nslice+1): # it is the NE layer when i = 1  
        #print(i)
        toclean = list(shapepoints[:,i,:]) # to remove repeating sublists
        cleaned = [list(i) for i in set(map(tuple, toclean))]
        # calculate volume of iregular made of n points in a 3D space
        with open(f'03_slice\\{seq}_slice_shell{i}_t{timeslice}.xvg', 'w') as outfile:
            np.savetxt(outfile, cleaned, fmt='%.4f')

    # delete repeating points
    print("--- %s seconds ---" % (time.time() - start_time))

# %%

#%%
# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
from operator import sub
from numpy import zeros, ones, arange, asarray, concatenate
import time
import os, glob
import threading

def min_angle(ne, pmpoint, center):
    '''
    find the minimal angle between two vectors center-to-p1 and center-to-p2
    '''
    tmp = []
    for i in range(0,len(ne)):
        # if pmpoint[0]-2 <= ne[i][0] <= pmpoint[0]+2:
        v1 = list(map(sub, ne[i], center))
        v2 = list(map(sub, pmpoint, center))
        # print(v1, v2)
        cross = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)),4) # to obtain the angle by math.acos
        tmp.append(cross)
        if cross == 1:
            return ne[i] # zero degree
            break
    return ne[tmp.index(max(tmp))]

shells = 7 

start_time = time.time() # Measure time elapsed in Python

# def cut_docked(ne,pmpoint,center, edge = 100):
#      for i in range(0,len(ne)):
#         v1 = list(map(sub, ne[i], center))
#         v2 = list(map(sub, pmpoint, center))
#         if np.linalg.norm(v2)-edge/31.3 > np.linalg.norm(v1): # pixel
#              pmpoint_new = center + v2*((np.linalg.norm(v2)-edge/31.3) / np.linalg.norm(v1))
             
     
'''
Read data file
'''

# os.chdir("I:\\Bing\\fluorescence\\Results\\SIM\\")   

# os.getcwd() 
# file = glob.glob(f'*.npy')

# mask_t = np.load(file) # T, z, C(PM/NE), x, y

file_dir = 'I:\\Bing\\fluorescence\\3D\\00_index\\'
output_dir = 'I:\\Bing\\fluorescence\\3D\\01_findpair\\'

file_list=[]
with open('../file_list.txt','r') as f:
	for line in f:
		file_list.append(line.strip('\n'))

for file_name in file_list[24:]:
    print(file_name)
    # raw_mask = np.load('masks\\' + seq + '.npy') # T, z, C(PM/NE), x, y
    k1 = np.loadtxt(f'{file_dir}{file_name}_bin2_ne_index.xvg')
    k2 = np.loadtxt(f'{file_dir}{file_name}_bin2_pm_index.xvg')

    center = [float(k1[0][0]),float(k1[0][1]),float(k1[0][2])] # Get center, ne, and pm

    ne=[]
    for i in range(1, len(k1)):
        ne.append([float(k1[i][0]),float(k1[i][1]),float(k1[i][2])])
        
    pm=[]
    for i in range(0, len(k2)):
        pm.append([float(k2[i][0]),float(k2[i][1]),float(k2[i][2])])
    break

'''
find the minimum angles between the vector of center and ne, and the vector of center and pm
'''

for shell in range(shells):
# def find_pair(shell):
    pair_out=open(f'{output_dir}{file_name}_bin2_pair_ne-pm_shell{shell}.xvg', 'w') # output
    t1 = int(len(pm) * (shell) / shells)
    t2 = int(len(pm) * (shell+1) / shells)
    print(f'{len(pm)}, {t2}')

    for i in range(t1,t2):
        #print(i)
        #pair_out.write(i)
        pair_out.write(' '.join(map(str, min_angle(ne, pm[i], center) + pm[i])))
        pair_out.write("\n") 

# threadpool = []

# for shell in range(shells):
#     th = threading.Thread(target=find_pair, args=(shell,))
#     threadpool.append(th)
# for th in threadpool:
#     th.start()
# for th in threadpool:
#     threading.Thread.join(th)

# pair_out=open(f'{seq}_pair_ne-pm.xvg', 'w') # output
# t1 = 0
# t2 = int(len(pm)/7.0)
# print(len(pm))
# print(t2)
# for i in range(t1,t2):
#     #print(i)
#     #pair_out.write(i)
#     pair_out.write(' '.join(map(str, min_angle(ne, pm[i], center) + pm[i])))
#     pair_out.write("\n")
#     #print(i, ' '.join(map(str, min_angle(ne, pm[i], center) + pm[i])) \n, end=" ", file=pair_out)

    

print("--- %s seconds ---" % (time.time() - start_time))

# %%

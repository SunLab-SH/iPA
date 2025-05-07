#%%
# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
from operator import sub
from numpy import zeros, ones, arange, asarray, concatenate
import time
import os, glob, time

def min_angle(ne, pmpoint, center):
    '''
    find the minimal angle between two vectors center-to-p1 and center-to-p2
    '''
    tmp = []
    for i in range(0,len(ne)):
        v1 = list(map(sub, ne[i], center))
        v2 = list(map(sub, pmpoint, center))
        cross = round(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)),4) # to obtain the angle by math.acos
        tmp.append(cross)
        if cross == 1:
            return ne[i] # zero degree
            break
    return ne[tmp.index(max(tmp))]

start_time = time.time() # Measure time elapsed in Python

'''
Read data file
'''
#%%
os.chdir("I:\\Bing\\fluorescence\\Results\\SIM\\")   

os.getcwd() 
# file = glob.glob(f'*.npy')

# mask_t = np.load(file) # T, z, C(PM/NE), x, y
seq = 'cut_0min_3_5_SIM'
print(seq)
# raw_mask = np.load('masks\\' + seq + '.npy') # T, z, C(PM/NE), x, y

k1 = np.loadtxt(f'{seq}_ne_index.xvg')
k2 = np.loadtxt(f'{seq}_pm_index.xvg')

center = [float(k1[0][0]),float(k1[0][1]),float(k1[0][2])] # Get center, ne, and pm

ne=[]
for i in range(1, len(k1)):
    ne.append([float(k1[i][0]),float(k1[i][1]),float(k1[i][2])])
    
pm=[]
for i in range(0, len(k2)):
    pm.append([float(k2[i][0]),float(k2[i][1]),float(k2[i][2])])

'''
find the minimum angles between the vector of center and ne, and the vector of center and pm
'''

for shell in range(7):
    pair_out=open(f'{seq}_pair_ne-pm_shell{shell}.xvg', 'w') # output
    t1 = int(len(pm) * (shell) /7.0)
    t2 = int(len(pm) * (shell+1) /7.0)
    print(len(pm))
    print(t2)
    for i in range(t1,t2):
        #print(i)
        #pair_out.write(i)
        pair_out.write(' '.join(map(str, min_angle(ne, pm[i], center) + pm[i])))
        pair_out.write("\n") 
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

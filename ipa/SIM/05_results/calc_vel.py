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

#%%
seqs = ['2.8-30_P3S3-1', '2.8-30_P3S3-2',
        '2.8-30_P10S2-1', '2.8-30_P10S2-2', '2.8-30_P16', 
        '2.8-30_P17S1-1', '2.8-30_test', '16.7-5_P4', 
        '16.7-5_P7-1', '16.7-5_P7-2', '16.7-5_P10S3-1', 
        '16.7-5_P10S3-2', '16.7-5_P20', '16.7-30_P14S1-1',
        '16.7-30_P14S1-2', '16.7-30_P14S1-3', '16.7-30_P21',
        '16.7-30_P22', '16.7-30_P22-1'] # '2.8-30_P2', 

N_shell = 17

# cellnumber = '2.8-30_P3S3-1'
for cellnumber in seqs:
    start_time = time.time() # Measure time elapsed in Python

    '''
    Read data file
    '''

    # data_dir = 
    with open('../../05_rdf/build/'+ str(cellnumber) +'_t8_velocity_dist.xvg') as fp:
        k1 = fp.read().splitlines()
        # k1 = ''.join(k1).strip('\n')

    N_isg = len(k1)/N_shell
    print(len(k1), N_isg)
    
    '''
    calculate the rdf profile
    '''

    isg_in_shell = []
    for i in range(N_isg):
        th = k1.index(i+' ')
        th/len(k1)
        isg = N_isg - k1[int(N_isg*1*i) : int(N_isg*1*(i+1))].count('\n')
        isg_in_shell.append(isg)
    
    isgout=open('../../05_rdf/results/'+ str(cellnumber) + '_isgrdf_ne-pm.xvg', 'w') # output

    isg = np.zeros(16)
    #print(k1[1,0])

    print(f'isg total number: {isg_in_shell[-1]}')

    for i in range(1, N_shell):
        isg[i-1] = abs(isg_in_shell[i]-isg_in_shell[i-1])/isg_in_shell[-1]

    print(f'mean: {np.mean(isg)}')

    np.savetxt(isgout, isg, fmt='%.4f')
    isgout.close()
    

print("--- %s seconds ---" % (time.time() - start_time))


# %%

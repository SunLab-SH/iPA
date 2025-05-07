# coding="utf-8"
#%%

import os, sys
import re
from tabnanny import check
from importlib_metadata import method_cache

filedir = os.path.dirname(os.path.abspath(__file__))
os.chdir(filedir)
sys.path.append(os.path.pardir)
print(os.getcwd())

import numpy as np
import threading
from common.parser import arg
from common.preprocess import *
import glob 
from multiprocessing.dummy import Pool as ThreadPool
from common import preprocess, distances_generator, outplot
from scipy.spatial.distance import cdist


# %%

def obtain_dir():
    for maindir, subdir, file_name_list in os.walk(os.path.join(arg.root_dir, arg.data_root_dir), topdown=False):
        locate = np.array(subdir)
    subdir1 = [  os.path.join(arg.root_dir, arg.data_root_dir, subdir) for subdir in locate]
    # print(subdir1)
    subdir2 = []
    for dir in subdir1:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            condition = np.array(subdir)
        for singlecondition in condition:
            subdir2.append( os.path.join(dir, singlecondition) )
    # print(subdir2)
    subdir3 = []
    for dir in subdir2:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            replica = np.array(subdir)
        for singlereplica in replica:
            subdir3.append( os.path.join(dir, singlereplica) )
    # print(subdir3)


    return subdir3 

subdir3 = obtain_dir()
print(subdir3)
# %%

id = [ (path.split('\\')[-1]) for path in subdir3  ]
print(id)
newid = []
for singleid in id:
    if '\uf00d' in singleid:
        newid.append(singleid.replace('\uf00d', ''))
    else: newid.append(singleid)
id =newid
print(id)
idpd = pd.DataFrame(columns=['ID'], data=id)

print(idpd)
idpd.to_csv('./dataid.csv', index=False)
# %%

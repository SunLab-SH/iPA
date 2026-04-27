#%%

import os
import cv2
import sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--baseDir", help="directory in which all your network will live", default='models')
parser.add_argument("--name", help="name of your network", default='N2V3D')
parser.add_argument("--dataPath", help="The path to your data")
parser.add_argument("--fileName", help="name of your data file", default="*.tif")
parser.add_argument("--output", help="The path to your data to be saved", default='predictions.tif')
parser.add_argument("--dims", help="dimensions of your data", default='YX')
parser.add_argument("--tile", help="will cut your image [TILE] times in every dimension to make it fit GPU memory", default=1, type=int)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# We import all our dependencies.
from n2v.models import N2V
from n2v.internals.N2V_DataGenerator import N2V_DataGenerator
import numpy as np
from matplotlib import pyplot as plt
from tifffile import imread
from tifffile import imwrite
from glob import glob
from os.path import join
from matplotlib.image import imread, imsave
from PIL import Image
# A previously trained model is loaded by creating a new N2V-object without providing a 'config'.
model_name = args.name
# the base directory in which our model will live
basedir = args.baseDir
model = N2V(config=None, name=model_name, basedir=basedir)


#%%
tiles = (args.tile, args.tile)

if 'Z' in args.dims or 'C' in args.dims:
    tiles = (1, args.tile, args.tile)

if 'Z' in args.dims and 'C' in args.dims:
    tiles = (1, args.tile, args.tile, 1)

#''' raw
datagen = N2V_DataGenerator()
imgs = datagen.load_imgs_from_directory(directory = args.dataPath, dims=args.dims, filter=args.fileName)
# files=glob(join(args.dataPath, args.fileName))
files = glob(join(args.dataPath, args.fileName))
files.sort()
names=[]
for f in files:
    names.append(f)
#'''

#''' HCY

imgs = os.listdir(args.dataPath)

#'''

print('started...')

for i, img in enumerate(imgs):
    # print(f'\rProcess {i}/{len(imgs)}', end='')

    if len(imgs) == 1:
        filename=args.output+args.fileName       
    if len(imgs) > 1:
        # filename=args.output+names[i][len(args.dataPath):]
        filename = os.path.join(args.dataPath, img)
    # if '783_5_200' not in filename:
    #     continue

    ''' raw
    img_=img[0,...]
    if len(img_.shape)>len(args.dims):
        img_=img_[...,0]
    #'''

    #''' hcy
    img_ = cv2.imread(filename, 0)
    print(f'\n{img_.shape}')
    #'''


    pred = model.predict( img_, axes=args.dims, n_tiles=tiles)
    
    print(f'pred: {pred.shape}')
    filename = os.path.join(args.output, img)
    cv2.imwrite(filename, pred)
    print(f'Output1: {filename}')

    imsave(filename, np.clip(pred,0.0,1.0),cmap='gray')
    img = Image.open(filename).convert('L')
    out_name = filename[:-3]+'png'
    print(f'Output2: {out_name}')
    #img.save(out_name)
    


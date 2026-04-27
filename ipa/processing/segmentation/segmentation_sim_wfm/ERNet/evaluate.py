
import math
import os

import torch
import torch.nn as nn
import time 

import torch.optim as optim
import torchvision
from torch.autograd import Variable

from .models import *
from .datahandler import *

import matplotlib.pyplot as plt
import glob
import numpy as np
from PIL import Image

def remove_dataparallel_wrapper(state_dict):
	r"""Converts a DataParallel model to a normal one by removing the "module."
	wrapper in the module dictionary

	Args:
		state_dict: a torch.nn.DataParallel state dictionary
	"""
	from collections import OrderedDict

	new_state_dict = OrderedDict()
	for k, vl in state_dict.items():
		name = k[7:] # remove 'module.' of DataParallel
		new_state_dict[name] = vl

	return new_state_dict


def changeColour(I): # change colours (used to match WEKA output)
    Inew = np.zeros(I.shape + (3,)).astype('uint8')
    for rowidx in range(I.shape[0]):
        for colidx in range(I.shape[1]):
            if I[rowidx][colidx] == 0:
                Inew[rowidx][colidx] = [198,118,255]
            elif I[rowidx][colidx] == 127:
                Inew[rowidx][colidx] = [79,255,130]
            elif I[rowidx][colidx] == 255:
                Inew[rowidx][colidx] = [255,0,0]
    return Inew


def EvaluateModel(opt):

    try:
        os.makedirs(opt.out)
    except IOError:
        pass

    opt.fid = open(opt.out + '/log.txt','w')
    print(opt)
    print(opt,'\n',file=opt.fid)
    
    net = GetModel(opt)

    if not opt.weights or not os.path.exists(opt.weights):
        print(f"Error: weights file '{opt.weights}' not found or not specified")
        print("Please specify a valid weights file path using --weights parameter")
        return
    
    checkpoint = torch.load(opt.weights)
    if opt.cpu:
        net.cpu()
    
    print('loading checkpoint',opt.weights)
    if opt.undomulti:
        checkpoint['state_dict'] = remove_dataparallel_wrapper(checkpoint['state_dict'])
    net.load_state_dict(checkpoint['state_dict'])

    if opt.root.split('.')[-1] == 'png' or opt.root.split('.')[-1] == 'jpg':
        imgs = [opt.root]
    else:
        imgs = []
        imgs.extend(glob.glob(opt.root + '/*.jpg'))
        imgs.extend(glob.glob(opt.root + '/*.png'))
        imgs.extend(glob.glob(opt.root + '/*.tif'))
        if len(imgs) == 0: # scan everything
            imgs.extend(glob.glob(opt.root + '/**/*.jpg',recursive=True))
            imgs.extend(glob.glob(opt.root + '/**/*.png',recursive=True))
            imgs.extend(glob.glob(opt.root + '/**/*.tif',recursive=True))

    imageSize = opt.imageSize

    for i, imgfile in enumerate(imgs):
        img = np.array(Image.open(imgfile))/255


        if len(img.shape) > 2:
            print('removing colour channel')
            img = img[:,:,0]

        print(np.min(img),np.max(img),img.shape)

        h,w = img.shape[0], img.shape[1]
        if imageSize == 0:
            imageSize = 250
            while imageSize+250 < h and imageSize+250 < w:
                imageSize += 250
            print('Set imageSize to',imageSize)


        # img_norm = (img - np.min(img)) / (np.max(img) - np.min(img)) 
        images = []

        images.append(img[:imageSize,:imageSize])
        images.append(img[h-imageSize:,:imageSize])
        images.append(img[:imageSize,w-imageSize:])
        images.append(img[h-imageSize:,w-imageSize:])

        proc_images = []
        for idx,sub_img in enumerate(images):
            pil_sub_img = Image.fromarray((sub_img*255).astype('uint8'))
            sub_tensor = toTensor(pil_sub_img)
            print('\r[%d/%d][%d/%d], shape is %dx%d - ' % (idx+1,len(images),i+1,len(imgs),sub_tensor.shape[1],sub_tensor.shape[2]),end='')
            sub_tensor = sub_tensor.unsqueeze(0)

            with torch.no_grad():
                if opt.cpu:
                    sr = net(sub_tensor)
                else:
                    sr = net(sub_tensor.cuda())
                sr = sr.cpu()

                m = nn.LogSoftmax(dim=0)
                sr = m(sr[0])
                sr = sr.argmax(dim=0, keepdim=True)
                

                pil_sr_img = toPIL(sr.float() / (opt.nch_out - 1))
                # pil_sr_img.save(opt.out + '/segmeneted_output_' + str(i) + '_' + str(idx) + '.png')
                # pil_sub_img.save(opt.out + '/imageinput_' + str(i) + '_' + str(idx) + '.png')

                proc_images.append(pil_sr_img)
            
        # stitch together
        img1 = proc_images[0]
        img2 = proc_images[1]
        img3 = proc_images[2]
        img4 = proc_images[3]

        woffset = (2*imageSize-w) // 2
        hoffset = (2*imageSize-h) // 2

        img1 = np.array(img1)[:imageSize-hoffset,:imageSize-woffset]
        img3 = np.array(img3)[:imageSize-hoffset,woffset:]
        top = np.concatenate((img1,img3),axis=1)

        img2 = np.array(img2)[hoffset:,:imageSize-woffset]
        img4 = np.array(img4)[hoffset:,woffset:]
        bot = np.concatenate((img2,img4),axis=1)

        oimg = np.concatenate((top,bot),axis=0)
        
        oimg[:10,:] = 0
        oimg[-10:,:] = 0
        oimg[:,:10] = 0
        oimg[:,-10:] = 0

        # oimg = changeColour(oimg)  // whether to use colours similar to the WEKA plugin (purple background)

        print(imgfile,i)

        if opt.out == 'root': # save next to orignal
            ext = imgfile.split('.')[-1]
            Image.fromarray(oimg).save(imgfile.replace('.' + ext,'_out.png'))
        else:
            ext = imgfile.split('.')[-1]
            filename = os.path.basename(imgfile).replace('.' + ext,'_out.png')
            Image.fromarray(oimg).save('%s/%s' % (opt.out,filename))


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Evaluate ERNet model')
    # Get script directory for relative paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    parser.add_argument('--root', type=str, default=os.path.join(script_dir, 'testdata/input'), help='Path to input images')
    # Look for common model file extensions in pretrained_model directory
    pretrained_dir = os.path.join(script_dir, 'pretrained_model')
    default_weights = pretrained_dir
    if os.path.isdir(pretrained_dir):
        # Look for .pth, .pt, or .pkl files
        for ext in ['*.pth', '*.pt', '*.pkl']:
            weight_files = glob.glob(os.path.join(pretrained_dir, ext))
            if weight_files:
                default_weights = weight_files[0]  # Use the first found
                break
    
    parser.add_argument('--weights', type=str, default=default_weights, help='Path to model weights')
    parser.add_argument('--imageSize', type=int, default=1000, help='Image size for processing')
    parser.add_argument('--out', type=str, default=os.path.join(script_dir, 'testdata/output'), help='Output directory')
    parser.add_argument('--model', type=str, default='rcan', help='Model type')
    parser.add_argument('--nch_in', type=int, default=1, help='Number of input channels')
    parser.add_argument('--nch_out', type=int, default=2, help='Number of output channels')
    parser.add_argument('--n_resgroups', type=int, default=5, help='Number of residual groups')
    parser.add_argument('--n_resblocks', type=int, default=10, help='Number of residual blocks')
    parser.add_argument('--cpu', action='store_true', help='Use CPU instead of GPU')
    parser.add_argument('--undomulti', action='store_true', help='Remove DataParallel wrapper')
    parser.add_argument('--n_feats', type=int, default=64, help='Number of features')
    parser.add_argument('--reduction', type=int, default=16, help='Reduction factor')
    parser.add_argument('--narch', type=int, default=0, help='Network architecture number')
    parser.add_argument('--multigpu', action='store_true', help='Use multiple GPUs')
    parser.add_argument('--dataset', type=str, default='imagedataset', help='Dataset type')
    parser.add_argument('--batchSize', type=int, default=16, help='Batch size')
    parser.add_argument('--batchSize_test', type=int, default=1, help='Test batch size')
    parser.add_argument('--lr', type=float, default=0.0001, help='Learning rate')
    parser.add_argument('--nepoch', type=int, default=10, help='Number of epochs')
    parser.add_argument('--ntrain', type=int, default=0, help='Number of training samples')
    parser.add_argument('--ntest', type=int, default=10, help='Number of test samples')
    parser.add_argument('--test', action='store_true', help='Test mode')
    parser.add_argument('--log', action='store_true', help='Enable logging')
    parser.add_argument('--saveinterval', type=int, default=10, help='Save interval')
    parser.add_argument('--testinterval', type=int, default=1, help='Test interval')
    parser.add_argument('--plotinterval', type=int, default=1, help='Plot interval')
    parser.add_argument('--scheduler', type=str, default='', help='Learning rate scheduler')
    parser.add_argument('--workers', type=int, default=6, help='Number of workers')
    
    opt = parser.parse_args()
    
    EvaluateModel(opt)






            




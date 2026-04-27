#%%

import argparse
import logging
import os

import numpy as np
import torch
import torch.nn.functional as F
from PIL import Image
from torchvision import transforms

from utils.data_loading import BasicDataset
from unet import UNet
from utils.utils import plot_img_and_mask
import mrcfile
import tifffile
import pandas as pd
import matplotlib.pyplot as plt


import os, sys
curpath = os.path.dirname(os.path.abspath(__file__))
os.chdir(curpath)
# os.path.pardir

sys.path.append(os.getcwd())
parpath = os.getcwd()
print(parpath)



#%%
def predict_img(net,
                full_img,
                device,
                scale_factor=1,
                out_threshold=0.5):
    net.eval()
    img = torch.from_numpy(BasicDataset.preprocess(full_img, scale_factor, is_mask=False))
    img = img.unsqueeze(0)
    img = img.to(device=device, dtype=torch.float32)

    with torch.no_grad():
        output = net(img)

        if net.n_classes > 1:
            probs = F.softmax(output, dim=1)[0]
        else:
            probs = torch.sigmoid(output)[0]

        tf = transforms.Compose([
            transforms.ToPILImage(),
            transforms.Resize((full_img.size[1], full_img.size[0])),
            transforms.ToTensor()
        ])

        full_mask = tf(probs.cpu()).squeeze()

    if net.n_classes == 1:
        return (full_mask > out_threshold).numpy()
    else:
        return F.one_hot(full_mask.argmax(dim=0), net.n_classes).permute(2, 0, 1).numpy()


def get_args():
    parser = argparse.ArgumentParser(description='Predict masks from input images')
    parser.add_argument('--model', '-m', default='vesicle_mask0.63_processed.pth', metavar='FILE',
                        help='Specify the file in which the model is stored')
    # parser.add_argument('--input', '-i', metavar='INPUT', nargs='+', help='Filenames of input images', required=True)
    # parser.add_argument('--output', '-o', metavar='OUTPUT', nargs='+', help='Filenames of output images')
    # parser.add_argument('--viz', '-v', action='store_true',
    #                     help='Visualize the images as they are processed')
    # parser.add_argument('--no-save', '-n', action='store_true', help='Do not save the output masks')
    parser.add_argument('--mask-threshold', '-t', type=float, default=0.5,
                        help='Minimum probability value to consider a mask pixel white')
    parser.add_argument('--scale', '-s', type=float, default=0.5,
                        help='Scale factor for the input images')
    parser.add_argument('--bilinear', action='store_true', default=False, help='Use bilinear upsampling')

    return parser.parse_args(args = [])


def get_output_filenames(args):
    def _generate_name(fn):
        return f'{os.path.splitext(fn)[0]}_OUT.png'

    return args.output or list(map(_generate_name, args.input))


def mask_to_image(mask: np.ndarray):
    if mask.ndim == 2:
        return Image.fromarray((mask * 255).astype(np.uint8))
    elif mask.ndim == 3:
        return Image.fromarray((np.argmax(mask, axis=0) * 255 / mask.shape[0]).astype(np.uint8))



# -------

def read_dataid():
    path = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\model_for_insulin_vesicle\\Pytorch-UNet\\utils\\dataid2_seg.txt'
    df = pd.read_csv(path, header=None)
    dataidlst = list(df.values)

    return [str(dataid[0]) for dataid in dataidlst]

# idlst = read_dataid()
# print(idlst)


def convert_img(img_3d):
    '3d normalization for 3d img'
    min_lac = np.min(img_3d)
    max_lac = 1 # np.max(img_3d)
    img_3d_new = (img_3d - min_lac )/ max_lac * 255


    # img_3d_new = np.zeros_like(img_3d)
    # for i in range(img_3d.shape[0]):
    #     curslice = img_3d[i]
    #     curslice = (curslice - np.min(curslice) )/ np.max(curslice) * 255
    #     img_3d_new[i] = curslice

    return img_3d_new.astype(np.uint8)

def get_filelist(dir, Filelist):
    newDir = dir
    if os.path.isfile(dir):
        Filelist.append(dir)
        # # 若只是要返回文件文，使用这个
        # Filelist.append(os.path.basename(dir))
    elif os.path.isdir(dir):
        for s in os.listdir(dir):
            # 如果需要忽略某些文件夹，使用以下代码
            #if s == "xxx":
                #continue
            newDir=os.path.join(dir,s)
            get_filelist(newDir, Filelist)
    return Filelist

#%%
def load_lac(idx):

    lacfilpath = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\Data\\raw_img\\Autosegmentation LIst.xlsx'
    lacc = 33.33
    df = pd.read_excel(lacfilpath)
    namelist = list(df['Cell ID'].values)
    laclst = list(df['LAC Normalizing Factor'].values)
    # print(namelist)

    for ii, name in enumerate(namelist):
            # print(idx, name, idx in name )
        if idx in name:
            lacc = laclst[ii]
    
    return lacc


def load_img(idx,):
    '''idx dataid'''
    # convert image
    imgpath = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\Data\\raw_img'

    namelist = get_filelist(imgpath, [])
    # print(namelist)
    for name in namelist:
        if f'{idx}' in name and 'pre_rec.tif' in name:
            # print( 'data name', name)
            img = tifffile.imread(name)
            break
        elif f'{idx}' in name and 'pre_rec.mrc' in name:
            # print( 'data name', name)
            img = mrcfile.open(name, permissive=True).data
            break
    
    lac = load_lac(idx)
    print(lac)
    img = img * lac
    # imgnew = copy.deepcopy(img)
    img[np.where(img > 1)] = 1

    print(np.min(img), np.max(img))
    return convert_img(img)


# 
# aa = load_img('TAK_30_2213_11-12')
# print(aa.shape )



#%%


def main():
    args = get_args()
    # in_files = args.input
    # out_files = get_output_filenames(args)

    dataidlst = read_dataid()
    dataidlst = ['islets_8_3-4', 'islets_9_2-3', 'islets_12_4']
    print(dataidlst)
    # dataidlst = dataidlst[:1]


    for idx in dataidlst:
        print(idx)
        net = UNet(n_channels=3, n_classes=2, bilinear=args.bilinear)

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logging.info(f'Loading model {args.model}')
        logging.info(f'Using device {device}')

        net.to(device=device)
        net.load_state_dict(torch.load(args.model, map_location=device))

        logging.info('Model loaded!')


        img3d = load_img(idx)
        imgshape = img3d.shape
        print(imgshape)
        outputmask = np.zeros_like(img3d)

        logging.info(f'\nPredicting image {idx} ...')

        for ii in range(imgshape[0]):

            img = img3d[ii]
            img_gray = np.zeros([img.shape[0],img.shape[1],3])
            img_gray[:,:,0], img_gray[:,:,1], img_gray[:,:,2] = img, img, img
            img_gray = img_gray.astype(np.uint8)
            # print(img.shape, img_gray.shape)
            img = Image.fromarray(img_gray)

            mask = predict_img(net=net,
                            full_img=img,
                            scale_factor=args.scale,
                            out_threshold=args.mask_threshold,
                            device=device)
            if ii % 100 == 0:
                print(ii, mask.shape, np.max(mask), np.min(mask))

            outputmask[ii] = mask[0] * 0 + mask[1] * 1
            # aa = Image.fromarray((mask * 255).astype(np.uint8))

        print(np.max(outputmask))
        outputfilename = f'F:\\salilab\\projects\\auto-segmentation_Aneesh\\model_for_insulin_vesicle\\Pytorch-UNet\\results\\isg_semantic_mask\\{idx}_pred_vesicle.tiff'
        tifffile.imsave(outputfilename, outputmask)
        print(f'Done with {idx}.')


if __name__ == '__main__':
    main()

# %%

import os
import sys
import argparse
import numpy as np
import tifffile
import skimage
from skimage.morphology import watershed
import json as j
from scipy.ndimage import binary_fill_holes
import pycocotools.mask as mask_util
from itertools import groupby
from schedule import Schedule


def init():
    parser = argparse.ArgumentParser()
    parser.add_argument('-tiffs', dest='tiffs', default='./tiff', type=str, help='directory of tiff')
    cfg = parser.parse_args()
    return cfg


def check(cfg):
    if not os.path.isdir(cfg.tiffs):
        print(f'Tiff directory not found! -> {cfg.tiffs}')
    print(f'Tiff directory: {cfg.tiffs}')
    if len(os.listdir(cfg.tiffs)) == 0:
        print(f'Directory is empty...Exit now!')
        exit()
    print('-------------------')


def label_mask(tif, flip):
    '''
    Author: Shuailin
    granule: (65280, 65280, 0)
    '''
    data = tif.astype(np.uint16)
    # data = np.flipud(data)    # flip?
    num, h, w, c = data.shape
    mask_label = np.zeros((num, h, w), dtype=np.uint8)
    
    assert c == 3
    data_flat = data.reshape(num*h*w, c)
    mask_label_flat = np.sum((data_flat == np.array([65280, 65280, 0])).astype(np.uint8), axis=-1)//3    
    mask_label = mask_label_flat.reshape(num, h, w).astype(np.uint8)
    
    if flip == 1:
        for i in range(num):
            mask_label[i] = np.flipud(mask_label[i])
            print(f'\rFlip Process: {i+1}/{num}', end='')
    
    return mask_label



def data_preprocess(data):
    '''
    Too slow and too stupid
    '''
    # uint16 -> uint8
    data = data.astype(np.uint16)
    # data = (data / 65535.)
    data = np.flipud(data)

    num, h,w,_ = data.shape
    new_data = np.zeros((num,h,w), dtype='i')
    for k in range(num):
        for i in range(h):
            for j in range(w):
                now = (round(data[k,i,j,0],3),round(data[k,i,j,1],3),round(data[k,i,j,2],3))
                if now == (45824, 52224, 52224): # background
                    new_data[k,i,j] = 0
                elif now == (8192, 32768, 48896): # kernel
                    new_data[k,i,j] = 0
                elif now == (57088, 57088, 57088): # membrane
                    new_data[k,i,j] = 0
                elif now == (57088, 32768, 48896): # mito
                    new_data[k,i,j] = 0
                elif now == (65280, 65280, 0): # granule
                    new_data[k,i,j] = 1
                else:
                    print(f'Data preprocess error: unknown raw data -> {now}')
                    exit(1)
                print(f'\rProcess: {k+1}/{num}  Row: {i+1}/{h}  Column: {j+1}/{w}', end='')
    return new_data

'''
backup
def get_bboxes(granule):
    granule,num = skimage.morphology.label(granule, return_num=True)

    if num == 0:
        return []

    bboxs = []
    for i in range(1, num+1):
        ins = np.zeros_like(granule)
        ins[(np.where(granule == i))] = 1
        ys, xs = np.where(granule == i)
        points = []
        for i in range(len(xs)):
            points.append((xs[i], ys[i]))
        left, top = sys.maxsize, sys.maxsize
        right, down = -1, -1
        for p in points:
            left = min(p[0], left)
            top = min(p[1], top)
            right = max(p[0], right)
            down = max(p[1], down)
        left -= 1
        top -= 1
        right += 1
        down += 1
        bbox = [int(left), int(top), int(right-left), int(down-top)]
        bboxs.append(bbox)
    return bboxs
'''

def is_around(mask, x, y):
    '''
      1 
    1 . 1
      1 
    '''
    h, w = mask.shape
    h -= 1
    w -= 1
    # top
    if x != 0 and mask[x-1][y] != 1:
        return False
    # down
    if x != h and mask[x+1][y] != 1:
        return False
    # left
    if y != 0 and mask[x][y-1] != 1:
        return False
    # right
    if y != w and mask[x][y+1] != 1:
        return False
    return True


def is_isolated(mask, x, y):
    '''
      .  
    . 1 .
      . 
    '''
    h, w = mask.shape
    h -= 1
    w -= 1
    # top
    if x != 0 and mask[x-1][y] == 1:
        return False
    # down
    if x != h and mask[x+1][y] == 1:
        return False
    # left
    if y != 0 and mask[x][y-1] == 1:
        return False
    # right
    if y != w and mask[x][y+1] == 1:
        return False
    return True


def is_clamped(mask, x, y):
    ''' 
    description: the point mask[x][y] must equal 1 and be isolated
    return: [top, right, down left]

    1 . 1
    . 1 .
    . . .

    '''

    h, w = mask.shape
    h -= 1
    w -= 1
    res = [False]*4
    # top
    if x != 0:
        if y != 0 and y != w:
            if mask[x-1][y-1] == 1 and mask[x-1][y+1] == 1:
                res[0] = True
    # down
    if x != h:
        if y != 0 and y != w:
            if mask[x+1][y-1] == 1 and mask[x+1][y+1] == 1:
                res[2] = True
    # left
    if y != 0:
        if x != 0 and x != h:
            if mask[x-1][y-1] == 1 and mask[x+1][y-1] == 1:
                res[3] = True
    # right
    if y != w:
        if x != 0 and x != h:
            if mask[x-1][y+1] == 1 and mask[x+1][y+1] == 1:
                res[1] = True
    return res


def is_connected(mask, x, y):
    '''
    description: the point mask[x][y] must equal 1 and be isolated and not be clamped
    return: [top right down left]

    . . .
    . 1 .
    1 . .

    '''
    h, w = mask.shape
    h -= 1
    w -= 1
    res = [False]*4
    # top-left
    if x != 0 and y != 0:
        if mask[x-1][y-1] == 1:
            res[0] = True
            res[3] = True
    # top-right
    if x != 0 and y != w:
        if mask[x-1][y+1] == 1:
            res[0] = True
            res[1] = True
    # down-right
    if x != h and y != w:
        if mask[x+1][y+1] == 1:
            res[1] = True
            res[2] = True
    # down-left
    if x != h and y != 0:
        if mask[x+1][y-1] == 1:
            res[2] = True
            res[3] = True
    return res



def fix_binary_mask(mask):
    res = mask.copy()
    res = binary_fill_holes(res).astype(int)
    xs,ys = np.where(mask == 1)
    assert len(xs) == len(ys)
    if len(xs) == 1:
        return res
    for i in range(len(xs)):
        x = xs[i]
        y = ys[i]
        if not is_isolated(mask, x, y):
            continue
        clamped = is_clamped(mask, x, y)
        if any(clamped):
            if clamped[0]:
                res[x-1][y] = 1
            if clamped[1]:
                res[x][y+1] = 1
            if clamped[2]:
                res[x+1][y] = 1
            if clamped[3]:
                res[x][y-1] = 1
            continue
        connected = is_connected(mask, x, y)
        try:
            assert any(connected) == True
        except AssertionError as e:
            print(f'\n\nlen xs = {len(xs)}')
            exit()
        if connected[0]:
            res[x-1][y] = 1
        if connected[1]:
            res[x][y+1] = 1
        if connected[2]:
            res[x+1][y] = 1
        if connected[3]:
            res[x][y-1] = 1
    res = binary_fill_holes(res).astype(int)
    return res


def binary_mask_to_rle(binary_mask):
    rle = {'counts': [], 'size': list(binary_mask.shape)}
    counts = rle.get('counts')
    for i, (value, elements) in enumerate(groupby(binary_mask.ravel(order='F'))):
        if i == 0 and value == 1:
            counts.append(0)
        counts.append(len(list(elements)))
    return rle


def get_bboxes_segs(granule):
    granule,num = skimage.morphology.label(granule, return_num=True)

    if num == 0:
        return [], []

    bboxs = []
    segs = []
    for i in range(1, num+1):
        ins = np.zeros_like(granule)
        ins[(np.where(granule == i))] = 1
        raw_ins = ins.copy()
        # generate segmentation mask (RLE)
        raw_ins = raw_ins.astype('uint8')
        now_seg = binary_mask_to_rle(raw_ins)
        segs.append(now_seg)

        # generate bbox (xywh)
        ys, xs = np.where(granule == i)
        points = []
        for i in range(len(xs)):
            points.append((xs[i], ys[i]))
        left, top = sys.maxsize, sys.maxsize
        right, down = -1, -1
        for p in points:
            left = min(p[0], left)
            top = min(p[1], top)
            right = max(p[0], right)
            down = max(p[1], down)
        left -= 1
        top -= 1
        right += 1
        down += 1
        bbox = [int(left), int(top), int(right-left+1), int(down-top+1)]
        bboxs.append(bbox)
    return bboxs, segs
    


def process_tiff(tiff, ins_label, flip):
    data = tifffile.imread(tiff)
    tiff = tiff.replace('_pre_rec_labels', '')
    data = label_mask(data, flip)
    print(f'Preprocess accomplished!\n', end='')
    num, h, w = data.shape
    schedule = Schedule(num)
    for k in range(num):
        temp = data.copy()
        temp = temp[k,:,:] # temp is one image
        bboxs, segs = get_bboxes_segs(temp)
        if len(bboxs) != 0:
            now_label = {
                'image_name': tiff.split('/')[-1].split('.')[0]+'_'+str(k)+'.png',
                'bbox': bboxs,
                'height': h,
                'width': w,
                'segmentation': segs
            }
            ins_label.append(now_label)
        # print(f'\rbbox Process: {k+1}/{num}', end='')
        schedule.watch()
    return data, ins_label


def process_tiffs(tiff_dir):
    tiffs = os.listdir(tiff_dir)
    ins_label = []
    fliptable = {
        '766_10' : 0,
        '766_7'  : 0,
        '783_12' : 1,
        '769_5'  : 1,
        '766_2'  : 1,
        '766_11' : 0,
        '769_7'  : 1,
        '842_13' : 1,
        '931_9' : 1,
        '931_14' : 1,
        '822_7'  : 1,
        '784_7'  : 1,
        '784_6'  : 1,
        '822_6'  : 1,
        '784_4'  : 1,
        '822_4'  : 1,
        '783_6'  : 1,
        '785_7'  : 0,

        '842_12' : 1,
        '766_5'  : 1,
        '783_5'  : 1,

        '766_8'  : 1,
        '842_17' : 1,
        '784_5'  : 1
    }
    for index, tiff in enumerate(tiffs, 1):
        print(f'Tiff process: {index}/{len(tiffs)}')
        flip = -1
        for cell in fliptable.keys():
            if cell in tiff:
                flip = fliptable[cell]
                break
        print(f'Now tiff: {tiff}  flip: {flip}')
        assert flip != -1
        if '.tif' in tiff:
            data, ins_label = process_tiff(os.path.join(tiff_dir, tiff), ins_label, flip)
            print('\n------------------------------\n',end='')
    # f = open('./json/ins_label.json', 'w', encoding='utf-8')
    f = open('/p300/ihuman/debug/temp1.json', 'w', encoding='utf-8')
    j.dump(ins_label, f)
    return ins_label

if __name__ == "__main__":
    cfg = init()
    check(cfg)
    process_tiffs(cfg)


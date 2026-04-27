import argparse
import cv2 as cv
import json as j
import pycocotools.mask as mask_util
from itertools import groupby

'''
Description: transform coco-format json file to crop coco-format json file. There will distinguish train or val or test because the crt-file will supply crop information.
'''

def init():
    parser = argparse.ArgumentParser()
    parser.add_argument('-raw', default='./json/nocrop_ins_label.json', help='file which will to be transformed')
    parser.add_argument('-crt', default='/group/xiangyi/Chuanyang/crop_info/test_info.json', help='assist information')
    parser.add_argument('-out', default='./json/crop_ins_label.json', help='name of save file')
    cfg = parser.parse_args()
    return cfg


def binary_mask_to_rle(binary_mask):
    rle = {'counts': [], 'size': list(binary_mask.shape)}
    counts = rle.get('counts')
    for i, (value, elements) in enumerate(groupby(binary_mask.ravel(order='F'))):
        if i == 0 and value == 1:
            counts.append(0)
        counts.append(len(list(elements)))
    return rle


def crop(data, crop_info):
    # raw_anno_file = cfg.raw
    criterion = crop_info
    out_anno_file = '/p300/ihuman/debug/temp3.json'

    # f = open(raw_anno_file, 'r')
    # data = j.load(f)
    f = open(criterion, 'r')
    c = j.load(f)

    print('LOAD-DATA ACCOMPLISHED!')

    # modify image width
    name2id = {}
    id2lr = {}
    for index, img in enumerate(data['images']):
        raw_img_name = img['file_name']
        name2id[raw_img_name] = img['id']
        for nimg in c['info']:
            if nimg['img_name'] == raw_img_name:
                assert nimg['raw_w'] == img['width']
                assert nimg['raw_h'] == img['height']
                new_w = nimg['raw_w']-nimg['left']-(nimg['raw_w']-nimg['right'])
                data['images'][index]['width'] = int(new_w)
                id2lr[img['id']] = (nimg['left'], nimg['right'])
                # w:5 left:1 right:3   0 [1 2] 3 4   5-1-(5-3) = 2

    # modify bbox data
    count = 0
    img_sum = len(c['info'])
    for img in c['info']:
        count += 1
        img_name = img['img_name']
        if img_name not in name2id:
            continue
        id = name2id[img_name]
        left,right = id2lr[id]
        for anno in data['annotations']:
            if anno['image_id'] == id:
                anno['bbox'][0] -= img['left']
                rle = anno['segmentation']
                compressed_rle = mask_util.frPyObjects(rle, rle.get('size')[0], rle.get('size')[1])
                bmask = mask_util.decode(compressed_rle)
                bmask = bmask[..., left:right]
                new_rle = binary_mask_to_rle(bmask)
                anno['segmentation'] = new_rle

        print('\rProcess: %d/%d  %.2f%%' % (count, img_sum, round(float(count) / img_sum * 100, 2)), end='')

    with open(out_anno_file, 'w') as f:
        j.dump(data, f)
    return data


if __name__ == "__main__":
    cfg = init()
    crop(cfg)



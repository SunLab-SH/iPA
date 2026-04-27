import argparse
import os
import json
import numpy as np
from schedule import Schedule


def hcy_rle2bmask(data):
    """
    :param data: {'rle': xxx, 'size': (h, w)} json
    :return: binary mask
    """
    rle_data = data['rle']
    size = data['size']
    matrix = []
    count = ''
    for char in rle_data:
        if char == 'a':
            matrix.extend([0]*int(count))
            count = ''
        elif char == 'b':
            matrix.extend([1] * int(count))
            count = ''
        else:
            count += char
    matrix = np.array(matrix).reshape(size)
    return matrix


parser = argparse.ArgumentParser(description='将IS Model结果还原为最原始大小映射反馈给ihuman做label修正')
parser.add_argument('-m', '--mode', default='', help='train or val or test')
args = parser.parse_args()

# 输入参数设定
args.mode = 'train'
crop_infopath = '/p300/ihuman/dataset/annotations/train/crop_info_train.json'
merge_labelpath = '/p300/ihuman/dataset/annotations/train/merge_train.json'

# 输出参数设定
out_dir = '/p300/ihuman/dataset/feedback'

assert args.mode != ''

if args.mode not in ['train', 'val', 'test']:
    raise ValueError

split_info = json.load(open(crop_infopath, 'r'))
merge_result = json.load(open(merge_labelpath, 'r'))

train_cells = [
    '766_10',
    '766_7' ,
    '783_12',
    '769_5' ,
    '766_2' ,
    '766_11',
    '769_7' ,
    '842_13',
    '931_9' ,
    '931_14',
    '822_7' ,
    '784_7' ,
    '784_6' ,
    '822_6' ,
    '784_4' ,
    '822_4' ,
    '783_6' ,
    '785_7'
]
val_cells = [
    '821_12',
    '766_5',
    '783_5'
]
test_cells = [
    '766_8',
    '842_17',
    '784_5'
]
cell2ids = {}
if args.mode == 'train':
    for item in train_cells:
        cell2ids[item] = []
if args.mode == 'val':
    for item in val_cells:
        cell2ids[item] = []
if args.mode == 'test':
    for item in test_cells:
        cell2ids[item] = []

# 构建id映射到图片名字和对应segmentations
id2name = {}
id2segs = {}
for img in merge_result['images']:
    id2name[img['id']] = img['file_name']
    id2segs[img['id']] = []
    find = False
    for cell_name in cell2ids.keys():
        if cell_name in img['file_name']:
            cell2ids[cell_name].append(img['id'])
            find = True
            break
    if find is False:
        raise ValueError
for anno in merge_result['annotations']:
    id2segs[anno['image_id']].append(anno['segmentation'])

# 构建name映射到切分信息
name2split_info = {}
for item in split_info['info']:
    name2split_info[item['img_name']] = item

for cell_name in cell2ids.keys():
    result_numpy_matrix = None
    result_order = []
    ids = cell2ids[cell_name]
    schedule = Schedule(ids, name=cell_name)
    for id in ids:
        name = id2name[id]
        segs = id2segs[id]
        info = name2split_info[name]
        assert len(segs) != 0
        mask = np.zeros(segs[0]['size'], dtype='uint8')
        instance_count = 1
        for seg in segs:
            decode_seg = hcy_rle2bmask(seg)
            assert decode_seg.shape == mask.shape
            mask[np.where(decode_seg == 1)] = instance_count
            instance_count += 1
        add_left = np.zeros((mask.shape[0], info['left']))
        add_right = np.zeros((mask.shape[0], info['raw_w'] - info['right']))
        mask = np.hstack((add_left, mask, add_right))
        if result_numpy_matrix is None:
            result_numpy_matrix = mask.reshape((1, mask.shape[0], -1)).copy()
        else:
            result_numpy_matrix = np.vstack((result_numpy_matrix, mask.reshape((1, mask.shape[0], -1)).copy()))
        result_order.append(name)
        schedule.watch()
    if not os.path.isdir(os.path.join(out_dir, cell_name)):
        os.system(f'mkdir -p {os.path.join(out_dir, cell_name)}')
    np.save(os.path.join(out_dir, cell_name, 'bmask.npy'), result_numpy_matrix)
    with open(os.path.join(out_dir, cell_name, 'order_list.json'), 'w') as f:
        json.dump(result_order, f)
    print('')


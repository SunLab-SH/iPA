# read npy and order_list, organized in order
# re-assign label to 3d ordered mask

import numpy as np
import json
import tifffile
import os
from traceback import print_exc
from schedule import Schedule

def generate_proposal(root_dir, dir_name):
    dir_path = root_dir + dir_name
    fuse_dir_path = dir_path
    # fuse_dir_path = '/root/PycharmProjects/1910_instance_segmentation/components/granule/output/'
    paths = {
        'npy': dir_path + '/bmask.npy',
        'order': dir_path + '/order_list.json',
        'fuse': fuse_dir_path + '/%s_3d.json' % dir_name,
        'save_npy': root_dir + 'output/%s_3d.npy' % dir_name,
        'save_tif': root_dir + 'output/%s_3d.tif' % dir_name,
    }

    mask = np.load(paths['npy']).astype(np.int)
    order_list = json.load(open(paths['order'], 'r'))
    fuse_list = json.load(open(paths['fuse']))
    order_ids = [int(it.split('.')[0].split('_')[-1]) for it in order_list]

    sorted_ids = sorted(order_ids)
    order_ids = [sorted_ids.index(it) for it in order_ids]  # real id order

    mask_clone = mask.copy()
    for ith, index in enumerate(order_ids):
        mask_clone[index] = mask[ith]
    mask_3d = np.zeros_like(mask_clone, dtype=np.uint16)
    schedule = Schedule(fuse_list, name=dir_name)
    for label_ith, _list in enumerate(fuse_list):
        for item in _list:
            slice_idx, label_idx = item.split('_')
            slice_idx, label_idx = sorted_ids.index(int(slice_idx)), int(label_idx)
            mask_3d[slice_idx][mask_clone[slice_idx] == label_idx] = label_ith + 1
        schedule.watch()

    # verify
    print('fused labels %d' % np.unique(mask_3d).max())
    # save
    # np.save(paths['save_npy'], mask_3d)
    tifffile.imsave(paths['save_tif'], mask_3d)
    return

def generate_order(root_dir, dir_name):
    dir_path = root_dir + dir_name
    order_path = dir_path + '/order_list.json'
    order_list = json.load(open(order_path, 'r'))
    order_ids = [int(it.split('.')[0].split('_')[-1]) for it in order_list]
    sorted_ids = sorted(order_ids)
    return sorted_ids

def main():
    root_dir = '/p300/ihuman/dataset/feedback-pure/fuse/FS_004/test/'
    print(f'root directory: {root_dir}')
    dir_names = sorted(os.listdir(root_dir))

    orders = {}

    for dir_name in dir_names:
        if ('_' in dir_name) and ('.' not in dir_name): #766_10
            try:
                generate_proposal(root_dir, dir_name)
                order = generate_order(root_dir, dir_name)
                orders[dir_name] = order
            except Exception as e:
                print_exc()
                print('%s failed' % dir_name)
                pass
                exit()
    dump_path = root_dir + 'output/order_dicts.json'
    json.dump(orders, open(dump_path, 'w'))


if __name__ == '__main__':
    main()

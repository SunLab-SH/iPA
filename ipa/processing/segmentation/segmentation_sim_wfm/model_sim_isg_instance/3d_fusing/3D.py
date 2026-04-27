import numpy as np
import json
from os.path import join
from os import listdir
from schedule import Schedule
from multiprocessing import Pool
from traceback import print_exc

root_dir = '/p300/ihuman/dataset/feedback-pure/fuse/FS_003/test'
cell_names = sorted(listdir(root_dir))
# cell_names = cell_names[12:18] #12, 18
print(f'root directory: {root_dir}')
print(cell_names)


threshold = 0.5


def mask_iou(mask_a, index_a, mask_b, index_b):
    global threshold
    smooth = 1e-5
    iou = float(((mask_a == index_a) & (mask_b == index_b)).sum()) / (((mask_a == index_a) | (mask_b == index_b)).sum() + smooth)
    if iou > threshold:
        return True
    else:
        return False


def mark(now_granule, granule_3d):
    for granule in granule_3d:
        if len(now_granule & granule) != 0:
            granule.update(now_granule)
            return granule_3d
    granule_3d.append(now_granule)
    return granule_3d


def to_3d(cell_name):
    granule_3d = []
    mask = np.load(join(root_dir, cell_name, 'bmask.npy'))
    name_order = json.load(open(join(root_dir, cell_name, 'order_list.json'), 'r'))
    order = {}
    for img_name in name_order:
        order[int(img_name.split('.')[0].split('_')[-1])] = (img_name, name_order.index(img_name))
    valid_index = list(order.keys())
    valid_index = sorted(valid_index)
    schedule = Schedule(len(valid_index)-2, name=cell_name)
    for index in range(1, len(valid_index) - 1):
        up_layer_index = valid_index[index - 1]
        middle_layer_index = valid_index[index]
        down_layer_index = valid_index[index + 1]

        up = mask[order[up_layer_index][1]]
        middle = mask[order[middle_layer_index][1]]
        down = mask[order[down_layer_index][1]]

        up_middle = []
        middle_down = []

        # 上中两两计算
        for up_index in range(1, int(np.max(up)) + 1):
            for middle_index in range(1, int(np.max(middle)) + 1):
                if mask_iou(up, up_index, middle, middle_index):
                    up_middle.append(set([f'{up_layer_index}_{up_index}', f'{middle_layer_index}_{middle_index}']))

        # 中下两两计算
        for middle_index in range(1, int(np.max(middle)) + 1):
            for down_index in range(1, int(np.max(down)) + 1):
                if mask_iou(down, down_index, middle, middle_index):
                    middle_down.append(set([f'{middle_layer_index}_{middle_index}', f'{down_layer_index}_{down_index}']))

        # 上中与中下匹配
        for ins_um in up_middle:
            for ins_md in middle_down:
                if len(ins_um & ins_md) == 1:
                    granule_3d = mark((ins_um | ins_md), granule_3d)

        schedule.watch()

    for i in range(len(granule_3d)):
        granule_3d[i] = list(granule_3d[i])

    with open(join(root_dir, cell_name, cell_name+'_3d.json'), 'w') as f:
        json.dump(granule_3d, f)
    print(f'granule 3d number: {len(granule_3d)}')


"""多进程版本"""
process_pool = Pool(processes=len(cell_names))
for item in cell_names:
    # to_3d(item)
    process_pool.apply_async(to_3d, (item,))
process_pool.close()
process_pool.join()


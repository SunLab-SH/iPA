import os
import json as j
import argparse
from utils import *
import traceback, sys

parser = argparse.ArgumentParser()
parser.add_argument('-mode', type=str, default='', help='handle test or train or val')
parser.add_argument('-info', type=str, default='/group/xiangyi/Chuanyang/ins_label_10_23/', help='directory to split information file')
parser.add_argument('-anno', type=str, default='/root/mask-rcnn/ins_label/train_json/fix_crop_ins_label_w-150_box-3.json', help='path to original annotation file before split')
parser.add_argument('-inferSeg', type=str, default='/root/upsnet/output/upsnet/coco/upsnet_resnet50_ihuman_4gpu_scratch/train-ihuman-sst/results/segmentations_coco_train-ihuman-sst_results.json')
parser.add_argument('-srcImg', type=str, default='/group/xiangyi/Chuanyang/ins_label_10_23/', help='小区域图片目录')
parser.add_argument('-save', type=str, default='/root/ihuman/merge/annotations/', help='Save directory of result')
parser.add_argument('-limitScore', type=float, default=0.5, help='Limitation to the segmentations')
parser.add_argument('-conflictIOU', dest='mergeScore',type=float, default=0.5, help='IoU between two mask greater than this parameter then will be merged')
args = parser.parse_args()

# 输入参数设定
args.mode = 'train'
args.info = '/p300/ihuman/dataset/annotations/train'
args.anno = '/p300/ihuman/dataset/annotations/train/crop_ins_label_train.json'
args.inferSeg = '/root/upsnet/output/upsnet/coco/upsnet_resnet50_ihuman_4gpu_scratch/train-ihuman-sst/results/segmentations_coco_train-ihuman-sst_results.json'
args.srcImg = '/p300/ihuman/dataset/split_images/train'
args.save = '/p300/ihuman/dataset/annotations/train/'
args.limitScore = 0.5
args.conflictIOU = 0.5

assert args.mode != ''

if args.mode == 'test' or args.mode == 'train' or args.mode == 'val':
    pass
else:
    assert 1+1 == 3

args.info = os.path.join(args.info, f'split_info_{args.mode}.json')
args.save = os.path.join(args.save, f'merge_{args.mode}.json')


def mask_iou(now_mask_location_f, exist_mask_f, conflict_index_f):
    xs_f, ys_f = np.where(exist_mask_f == conflict_index_f)
    exist_mask_location_f = []
    for item_f in zip(xs_f, ys_f):
        exist_mask_location_f.append(item_f)
    joint = float(len(now_mask_location_f.intersection(exist_mask_location_f)))
    result = (joint/len(now_mask_location_f), joint/len(exist_mask_location_f))
    return max(result), True if result[0] == max(result) else False

print('load annotations...', end='')
global_split_info = j.load(open(args.info, 'r'))
anno = j.load(open(args.anno, 'r'))
infer_segs = j.load(open(args.inferSeg, 'r'))
img_names = os.listdir(args.srcImg)
print('\rload successfully!')

# construct split id maps to it's own segmentation
split_id2segs = {}
for seg in infer_segs:
    if seg['image_id'] not in split_id2segs:
        split_id2segs[seg['image_id']] = [seg]
    else:
        split_id2segs[seg['image_id']].append(seg)

# construct original id maps to it's corresponding images which are split.
origin_id2imgs = {}
for img_name in img_names:
    temp = img_name.strip().split('.')[0].split('_')
    origin_id = int(temp[1])
    split_id = int(temp[2])
    if origin_id not in origin_id2imgs:
        origin_id2imgs[origin_id] = [{
            'img_name': img_name,
            'split_id': split_id,
            'segmentations': split_id2segs[split_id]
        }]
    else:
        origin_id2imgs[int(temp[1])].append({
            'img_name': img_name,
            'split_id': int(temp[2]),
            'segmentations': split_id2segs[split_id]
        })

# construct original id maps to it's own information entailing name and size
origin_id2info = {}
for img in anno['images']:
    origin_id2info[img['id']] = {
        'img_name': img['file_name'],
        'size': (img['height'], img['width'])
    }


# main process of merge
merge_anno = {
    'images': [],
    'annotations': [],
    'type': 'instance',
    'categories': [{'supercategory': 'none', 'id': 1, 'name': 'granule'}]
}
conflict_counter = 0
schedule = Schedule(origin_id2imgs, name='Image')
ci_count = 1
cs_count = 1
gi_count = 1
gs_count = 1
for index, origin_id in enumerate(origin_id2imgs, 1):
    step1 = False
    step2 = False
    step3 = False
    step4 = False
    try:
        origin_name = origin_id2info[origin_id]['img_name']
        origin_size = origin_id2info[origin_id]['size']
        corresponding_imgs = origin_id2imgs[origin_id]
        origin_binary_mask = np.zeros(origin_size, dtype='uint8')

        step1 = True
        # collect all segs in original image first
        corresponding_segs = []
        for img in corresponding_imgs:
            expand_x, expand_y = global_split_info[img['img_name']]
            segs = img['segmentations']
            for seg in segs:
                if seg['score'] < args.limitScore:
                    continue
                corresponding_segs.append({
                    'score': seg['score'],
                    'segmentation': seg['segmentation'],
                    'expand': (expand_x, expand_y)
                })

        step2 = True

        # 全部采用新排序方法：
        index_copy = []
        for corresponding_seg in corresponding_segs:
            index_copy.append((corresponding_seg['score'], corresponding_segs.index(corresponding_seg)))
        index_copy = sorted(index_copy, key=lambda x: (x[0], x[1]), reverse=True)
        corresponding_segs = [corresponding_segs[item[1]] for item in index_copy]

        # handle each segmentation in list sorted by score
        # corresponding_segs = sorted(corresponding_segs, key=lambda x: (x['score'], x['segmentation'], x['expand']), reverse=True)

        instance_number = len(corresponding_segs)
        for seg_index, now_seg in enumerate(corresponding_segs, 1):
            exist_mask_location = set()
            x_cords, y_cords = np.where(origin_binary_mask != 0)
            for item in zip(x_cords, y_cords):
                exist_mask_location.add(item)

            expand_x, expand_y = now_seg['expand']
            now_binary_mask = mask_util.decode(now_seg['segmentation'])

            now_mask_location = set()
            x_cords, y_cords = np.where(now_binary_mask == 1)
            for x, y in zip(x_cords, y_cords):
                now_mask_location.add((x+expand_x, y+expand_y))

            if len(exist_mask_location & now_mask_location) == 0:
                for x, y in now_mask_location:
                    origin_binary_mask[x][y] = seg_index
            else:
                conflict_index_lis = set()
                for item in (exist_mask_location & now_mask_location):
                    row, col = item
                    conflict_index_lis.add(origin_binary_mask[row][col])
                assert 0 not in conflict_index_lis
                conflict_index_lis = list(conflict_index_lis)
                conflict_counter += 1

                # 冲突解决策略
                if len(conflict_index_lis) == 1:
                    # 两个mask冲突
                    conflict_index = conflict_index_lis[0]
                    max_iou, order = mask_iou(now_mask_location, origin_binary_mask, conflict_index)
                    if max_iou > args.mergeScore:
                        # 较大值大于阈值，合并两个mask
                        if order:
                            # 当前mask重叠部分占比为较大值，舍弃当前mask
                            pass
                        else:
                            # 冲突mask重叠部分占比为较大值，舍弃冲突mask
                            origin_binary_mask[np.where(origin_binary_mask == conflict_index)] = 0
                            for x, y in now_mask_location:
                                origin_binary_mask[x][y] = seg_index
                    else:
                        # 较大值小于阈值，保留两个mask且当前mask的score一定较小
                        remain_cords = now_mask_location.difference(exist_mask_location)
                        for x, y in remain_cords:
                            origin_binary_mask[x][y] = seg_index
                else:
                    # 与多个mask同时冲突
                    max_iou_lis = []
                    order_lis = []
                    for conflict_index in conflict_index_lis:
                        max_iou, order = mask_iou(now_mask_location, origin_binary_mask, conflict_index)
                        max_iou_lis.append(max_iou)
                        order_lis.append((order, conflict_index))
                    greater_index_lis = [(max_iou_lis.index(i), i) for i in max_iou_lis if i > args.mergeScore]
                    if len(greater_index_lis) > 0:
                        greater_index_lis = sorted(greater_index_lis, key=lambda x: (x[1], x[0]), reverse=True)
                        # 存在有冲突mask使得可以合并，挨着合并直到列表为空或者当前mask在合并中被丢弃则停止
                        now_disappear = False
                        while len(greater_index_lis) > 0 and now_disappear is False:
                            if order_lis[greater_index_lis[0][0]][0]:
                                # 当前mask重叠部分占比为较大值，舍弃当前mask
                                now_disappear = True
                            else:
                                # 冲突mask重叠部分占比为较大值，舍弃冲突mask
                                origin_binary_mask[
                                    np.where(origin_binary_mask == order_lis[greater_index_lis[0][0]][1])] = 0
                                exist_mask_location = set()
                                x_cords, y_cords = np.where(origin_binary_mask != 0)
                                for item in zip(x_cords, y_cords):
                                    exist_mask_location.add(item)
                            greater_index_lis.remove(greater_index_lis[0])
                        if now_disappear:
                            pass
                        else:
                            remain_cords = now_mask_location.difference(exist_mask_location)
                            for x, y in remain_cords:
                                origin_binary_mask[x][y] = seg_index
                    else:
                        # 与所有冲突mask的iou都小于阈值，那么保留当前mask且移除所有冲突位置
                        remain_cords = now_mask_location.difference(exist_mask_location)
                        for x, y in remain_cords:
                            origin_binary_mask[x][y] = seg_index

        step3 = True

        # instance_mask, instance_number = label(origin_binary_mask, return_num=True)
        generate_img = {
            'file_name': origin_name,
            'id': origin_id,
            'height': origin_size[0],
            'width': origin_size[1],
            'instances': instance_number
        }
        generate_segs = []
        for i in range(1, instance_number+1):
            one_instance_mask = np.zeros_like(origin_binary_mask)
            one_instance_mask[np.where(origin_binary_mask == i)] = 1
            temp_xs, temp_ys = np.where(one_instance_mask == 1)
            if len(temp_xs) == 0 and len(temp_ys) == 0:
                continue
            compressed_mask = hcy_bmask2rle(one_instance_mask)
            generate_seg = {
                'image_id': origin_id,
                'segmentation': compressed_mask,
            }
            generate_segs.append(generate_seg)

        step4 = True

        merge_anno['images'].append(generate_img)
        merge_anno['annotations'].extend(generate_segs)
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        error = str(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        with open('/p300/ihuman/debug/merge.log', 'a+', encoding='utf-8') as f:
            f.write(f'origin-id: {str(origin_id)}  index: {str(index)}\n')
            if step1:
                f.write(f'step1 >>\norigin-name: {str(origin_name)}\norigin-size: {str(origin_size)}\ncorresponding-imgs: ci{str(ci_count)}.json\n')
                with open(f'/p300/ihuman/debug/ci{str(ci_count)}.json', 'w') as fff:
                    j.dump(corresponding_imgs, fff)
                ci_count += 1
            if step2:
                f.write(f'step2 >>\ncorresponding-segs: cs{str(cs_count)}.json\n')
                with open(f'/p300/ihuman/debug/cs{str(cs_count)}.json', 'w') as fff:
                    j.dump(corresponding_segs, fff)
                cs_count += 1
            if step3:
                f.write(f'step3 >> <<\n')
            if step4:
                f.write(f'step4 >>\ngenerate-img: gi{str(generate_img)}.json\ngenerate-segs: gs{str(gs_count)}.json\n')
                with open(f'/p300/ihuman/debug/gi{str(gi_count)}.json', 'w') as fff:
                    j.dump(generate_img, fff)
                gi_count += 1
                with open(f'/p300/ihuman/debug/gs{str(gs_count)}.json', 'w') as fff:
                    j.dump(generate_segs, fff)
                gs_count += 1
            f.write(error+'\n')
            f.write('-------------------------\n')
        with open('/p300/ihuman/debug/temp.json', 'w', encoding='utf-8') as f:
            j.dump(merge_anno, f)
    schedule.watch()
    
merge_anno['conflict_number'] = conflict_counter

with open(args.save, 'w', encoding='utf-8') as f:
    j.dump(merge_anno, f)


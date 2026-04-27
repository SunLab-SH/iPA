from os.path import join, isdir, isfile
from os import system, listdir
import json
import random
import cv2
from schedule import Schedule
import numpy as np
import pycocotools.mask as mask_util

COLORS = [
    (31, 119, 180), (255, 127, 14), (44, 160, 44), (148, 103, 189),
    (227, 119, 194), (188, 189, 34), (23, 190, 207)
]

def rle_decode(data):
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

"""
用于验证小区域划分的结果正确与否，输入为大图数据集、大图标注文件、小区域图片数据集和小区域标注文件，输出为一张大图可视化与其对应的所有小区域可视化
注-1：大图为cropped & denoised，denoise的时候标注文件无变动，因此大图标注文件采用crop结果标注文件
注-2：在split小区域的时候加入了bbox >= 3的筛选条件，因此存在部分大图经过split后实际上没有对应的小区域图片，因此选择可视化对象的时候要从小区域标注文件中选
注-3：默认bbox与segment均可视化
"""

# 可选参数
draw_bbox = True
draw_segmentation = True
label_type = 'poly' # "rle" or "poly"

# 输入参数设定
denoise_imgdir = '/p300/ihuman/dataset/denoise_images/train'
denoise_labelpath = f'/p300/ihuman/dataset/annotations/train/crop_ins_label_train.json'
split_imgdir = '/p300/ihuman/dataset/split_images/train'
split_labelpath = f'/p300/ihuman/dataset/annotations/train/split_ins_label_{label_type}_train.json'

# 输出参数设定
out_dir = '/p300/ihuman/debug/split'
if isdir(out_dir):
    if len(listdir(out_dir)) != 0:
        system(f'rm {join(out_dir, "*")}')
else:
    print('output directory not found!')
    exit()
if isfile(denoise_labelpath) and isfile(split_labelpath):
    pass
else:
    print('label file not found!')
    exit()
if label_type != 'rle' and label_type != 'poly':
    print('label_type must set to "rle" or "poly"! ')
    exit()

# 执行逻辑
denoise_label = json.load(open(denoise_labelpath, 'r'))
split_label = json.load(open(split_labelpath, 'r'))

# 从小区域标注文件中选取可视化对象的大图id
denoise_ids = set()
for img in split_label['images']:
    denoise_ids.add(int(img['image_name'].split('_')[0]))
# denoise_id = random.choice(list(denoise_ids))
denoise_id = 3820

# 可视化大图
denoise_img = None
denoise_imgname = ''
denoise_bboxes = []
denoise_segs = []
for img in denoise_label['images']:
    if img['id'] == denoise_id:
        denoise_imgname = img['file_name']
        denoise_img = cv2.imread(join(denoise_imgdir, denoise_imgname))
        assert denoise_img is not None
for label in denoise_label['annotations']:
    if label['image_id'] == denoise_id:
        denoise_bboxes.append(label['bbox'])
        denoise_segs.append(label['segmentation'])
assert len(denoise_bboxes) != 0
# 画bbox
if draw_bbox:
    for bbox in denoise_bboxes:
        x, y, w, h = bbox
        cv2.rectangle(denoise_img, (int(x), int(y)), (int(x + w - 1), int(y + h - 1)), (0, 0, 255), 1)
# 画segmentation
if draw_segmentation:
    for seg in denoise_segs:
        compressed = mask_util.frPyObjects(seg, seg['size'][0], seg['size'][1])
        bmask = mask_util.decode(compressed)
        points = np.where(bmask == 1)
        xs = points[0]
        ys = points[1]
        color = random.choice(COLORS)
        for x, y in zip(xs, ys):
            denoise_img[x][y][0] = color[2]
            denoise_img[x][y][1] = color[1]
            denoise_img[x][y][2] = color[0]
cv2.imwrite(join(out_dir, denoise_imgname), denoise_img)
print('Origin image visualize successfully!')

# 获取所有对应的小区域图片与其对应的所有bbox
split_id2name = {}
split_bboxes = {}
split_segs = {}
for img in split_label['images']:
    if int(img['image_name'].split('_')[0]) == denoise_id:
        split_id2name[img['id']] = img['file_name']
        split_bboxes[img['id']] = []
        split_segs[img['id']] = []
for label in split_label['annotations']:
    if label['image_id'] in split_bboxes:
        split_bboxes[label['image_id']].append(label['bbox'])
        split_segs[label['image_id']].append(label['segmentation'])

# 可视化小区域图片
schedule = Schedule(split_id2name)
for img_id in split_id2name.keys():
    split_imgname = split_id2name[img_id]
    split_img = cv2.imread(join(split_imgdir, split_imgname))
    assert split_img is not None
    # 画bbox
    if draw_bbox:
        for bbox in split_bboxes[img_id]:
            x, y, w, h = bbox
            cv2.rectangle(split_img, (int(x), int(y)), (int(x + w - 1), int(y + h - 1)), (0, 0, 255), 1)
    # 画segmentation
    if draw_segmentation:
        for seg in split_segs[img_id]:
            if label_type == 'rle':
                compressed = mask_util.frPyObjects(seg, seg['size'][0], seg['size'][1])
                bmask = mask_util.decode(compressed)
                points = np.where(bmask == 1)
                xs = points[0]
                ys = points[1]
                color = random.choice(COLORS)
                for x, y in zip(xs, ys):
                    split_img[x][y][0] = color[2]
                    split_img[x][y][1] = color[1]
                    split_img[x][y][2] = color[0]
            elif label_type == 'poly':
                compressed = mask_util.frPyObjects(seg, split_img.shape[0], split_img.shape[1])
                bmask = mask_util.decode(compressed)
                points = np.where(bmask == 1)
                xs = points[0]
                ys = points[1]
                color = random.choice(COLORS)
                for x, y in zip(xs, ys):
                    split_img[x][y][0] = color[2]
                    split_img[x][y][1] = color[1]
                    split_img[x][y][2] = color[0]
    cv2.imwrite(join(out_dir, split_imgname), split_img)
    schedule.watch()


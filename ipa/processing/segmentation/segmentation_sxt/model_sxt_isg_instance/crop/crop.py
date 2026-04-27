from os.path import join
from os import listdir
import cv2
import json
import sys
import os

"""
对原图进行crop操作，输入有原始图片与标注图片
输出cropped图片与crop information
"""

# 参数
output_crop_img = True  # 是否保存crop后的img
output_crop_info = True  # 是否保存crop的信息

print(f'output_crop_img: {output_crop_img}')
print(f'output_crop_info: {output_crop_info}')
print('')

# 设定 - 修改为本地路径
label_dir = r'd:\Gitspace\ipa_full\data\labels\train'  # 标注图像目录
raw_imgdir = r'd:\Gitspace\ipa_full\data\raw_images\train'  # 原始图像目录
crop_imgdir = r'd:\Gitspace\ipa_full\data\crop_images\train'  # (可选)中间输出
cell_pattern = ''  # 可以设置特定的cell pattern过滤
crop_outdir = r'd:\Gitspace\ipa_full\iPA\examples\test_isg model\instance\data\crop_images\train'  # 输出cropped图像
anno_outpath = r'd:\Gitspace\ipa_full\iPA\examples\test_isg model\instance\crop\crop_info_train.json'  # 输出crop信息


# 获取931_14所有相关数据file path
cells = []
for i in listdir(label_dir):
    if cell_pattern in i:
        cells.append(i)


# 执行逻辑
# 1-1 处理label
info = {
    'info': []
}

# 创建输出目录
os.makedirs(crop_outdir, exist_ok=True)
os.makedirs(os.path.dirname(anno_outpath), exist_ok=True)

print(f'Processing {len(cells)} files...')

for count, file_name in enumerate(cells, 1):
    print(f'Processing {count}/{len(cells)}: {file_name}')
    now_info = {
        'img_name': file_name
    }
    label_path = join(label_dir, file_name)
    la = cv2.imread(label_path, cv2.IMREAD_GRAYSCALE)
    
    if la is None:
        print(f'Warning: Could not read label image {label_path}')
        continue
        
    h, w = la.shape
    now_info['raw_h'] = h
    now_info['raw_w'] = w
    left = sys.maxsize
    right = -1
    for i in range(w):
        temp = la[:, i]
        if any(temp):
            left = i
            break
    for i in reversed(range(w)):
        temp = la[:, i]
        if any(temp):
            right = i+1
            break
    
    # crop image
    raw_imgpath = join(raw_imgdir, file_name)
    raw_img = cv2.imread(raw_imgpath, cv2.IMREAD_GRAYSCALE)
    
    if raw_img is None:
        print(f'Warning: Could not read raw image {raw_imgpath}')
        continue
        
    raw_img = raw_img[:, left:right]
    now_info['left'] = left
    now_info['right'] = right
    
    # save cropped image
    crop_imgpath = join(crop_outdir, file_name)
    if output_crop_img:
        cv2.imwrite(crop_imgpath, raw_img)
    info['info'].append(now_info)

if output_crop_info:
    with open(anno_outpath, 'w', encoding='utf-8') as f:
        json.dump(info, f)


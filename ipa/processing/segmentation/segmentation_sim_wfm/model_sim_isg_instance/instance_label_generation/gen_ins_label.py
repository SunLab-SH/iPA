import argparse
import sys

sys.path.append('/p300/ihuman/code/ins_label/integrate')
sys.path.append('/p300/ihuman/code/ins_label')

from gen_bbox import process_tiffs
from trans2coco import bbox2coco
from trans2crop import crop
import json
import random
from os import listdir, system
from os.path import join
import cv2
from schedule import Schedule

# 最小宽度筛选
min_width = 150
# 最少bbox筛选
min_bbox_num = 3

train_cell = [
    '766_10',
    '766_7',
    '783_12',
    '769_5',
    '766_2',
    '766_11',
    '769_7',
    '842_13',
    '931_9',
    '931_14',
    '822_7',
    '784_7',
    '784_6',
    '822_6',
    '784_4',
    '822_4',
    '783_6',
    '785_7',
]


def init():
    parser = argparse.ArgumentParser(description='输入semantic tiff，输出instance label (coco format)')
    # 生成label参数
    parser.add_argument('--tiffs', dest='tiff_dir', default='', type=str, help='输入sem tiff文件夹目录')
    parser.add_argument('--crop-info', dest='crop_info', default='', help='crop information')
    # 验证结果参数
    parser.add_argument('--img-dir', dest='img_dir', default='', help='cropped image目录')
    parser.add_argument('--out-dir', dest='out_dir', default='', help='验证结果输出目录')
    cfg = parser.parse_args()

    assert cfg.tiff_dir != ''
    assert cfg.crop_info != ''
    assert cfg.img_dir != ''
    assert cfg.out_dir != ''
    return cfg


def main(tiff_dir, crop_info, crop_imgdir):
    # 第一步：连通度算法生成instance level无格式的bbox list
    # Debug: temp1.json
    ins_label = process_tiffs(tiff_dir)

    # 第二步：bbox list转换成coco标注格式
    # Debug: temp2.json
    ins_label = bbox2coco(ins_label)

    # 第三步：根据原图crop信息将coco bbox映射到crop image上
    # Debug: temp3.json
    ins_label = crop(ins_label, crop_info)

    # 第四步：清理两边互相多余的部分完成标注与图片完全匹配
    ins_label = clean(ins_label, crop_imgdir)

    # 第五步：清理宽度小于150或gt-bbox数量小于3的图片与标注
    ins_label = filter(ins_label, crop_imgdir)

    with open('/p300/ihuman/dataset/annotations/train/crop_ins_label_train.json', 'w', encoding='utf-8') as f:
        json.dump(ins_label, f)

    return ins_label


def verify(anno, crop_imgdir, out_imgdir):
    # 每一组选择190-210之间的一个来验证
    special = random.randint(190, 210)
    id_patterns = []
    for cell_id in train_cell:
        id_patterns.append(cell_id + '_' + str(special))
    all_imgnames = listdir(crop_imgdir)
    # 遍历所有cropped images
    schedule = Schedule(all_imgnames)
    for img_name in all_imgnames:
        # 遍历所有还没验证的cell组特定名字串
        for id_pattern in id_patterns:
            # 监测 cell-id+special 是否在名字中
            if id_pattern in img_name:
                img = cv2.imread(join(crop_imgdir, img_name))
                # 确定id
                now_id = -1
                for item in anno['images']:
                    if item['file_name'] == img_name:
                        now_id = item['id']
                        break
                assert now_id != -1
                # 获取这张image的所有cropped bbox
                bboxes = []
                for item in anno['annotations']:
                    if item['image_id'] == now_id:
                        bboxes.append(item['bbox'])
                assert len(bboxes) != 0
                for bbox in bboxes:
                    x, y, w, h = bbox
                    cv2.rectangle(img, (int(x), int(y)), (int(x + w - 1), int(y + h - 1)), (0, 0, 255), 1)
                cv2.imwrite(join(out_imgdir, img_name), img)
        schedule.watch()


def clean(anno, crop_imgdir):
    all_imgnames = listdir(crop_imgdir)
    anno_imgnames = []
    for i in anno['images']:
        anno_imgnames.append(i['file_name'])
    # 先清理没有anno，有image
    for i in all_imgnames:
        if i not in anno_imgnames:
            system(f'rm {join(crop_imgdir, i)}')
    # 清理有anno，但没有image TODO:应该改为补充上这部分image
    delete_imgs = []
    delete_imgid = []
    for i in anno['images']:
        if i['file_name'] not in all_imgnames:
            delete_imgs.append(i)
            delete_imgid.append(i['id'])
    for i in delete_imgs:
        anno['images'].remove(i)
    delete_bboxes = []
    for i in anno['annotations']:
        if i['image_id'] in delete_imgid:
            delete_bboxes.append(i)
    for i in delete_bboxes:
        anno['annotations'].remove(i)

    return anno


def filter(anno, crop_imgdir):
    global min_width, min_bbox_num

    # 根据筛选条件获取要删除img的id列表
    delete_imgids = []
    # 宽度筛选
    for img in anno['images']:
        if img['width'] < min_width:
            delete_imgids.append(img['id'])
    # bbox数量筛选
    id2bbox = {}
    for bbox in anno['annotations']:
        if bbox['image_id'] not in id2bbox:
            id2bbox[bbox['image_id']] = [bbox]
        else:
            id2bbox[bbox['image_id']].append(bbox)
    for id in id2bbox.keys():
        if len(id2bbox[id]) < min_bbox_num:
            delete_imgids.append(id)

    # 根据id获取要删除的img
    delete_imgs = []
    for img in anno['images']:
        if img['id'] in delete_imgids:
            delete_imgs.append(img)

    # 删除图片
    schedule = Schedule(delete_imgs)
    for img in delete_imgs:
        system(f'rm {join(crop_imgdir, img["file_name"])}')
        schedule.watch()
    for item in delete_imgs:
        anno['images'].remove(item)

    # 根据id获取要删除的anno
    delete_bboxes = []
    for bbox in anno['annotations']:
        if bbox['image_id'] in delete_imgids:
            delete_bboxes.append(bbox)

    # 删除标注
    for item in delete_bboxes:
        anno['annotations'].remove(item)

    return anno


if __name__ == '__main__':
    cfg = init()
    anno = main(cfg.tiff_dir, cfg.crop_info, cfg.img_dir)
    verify(anno, cfg.img_dir, cfg.out_dir)


import os
import json as j

'''
Description: transform ins_label.json to coco-format json file. There will not distinguish train or val or test. Just excute it!
'''

def process():
    data = j.load(open('./train_json/ins_label.json', 'r'))
    ins_label = {
        'images': [],
        'type': 'instance',
        'annotations': [],
        'categories': [{'supercategory': 'none', 'id': 1, 'name': 'granule'}]
    }
    image_id = 1
    anno_id = 1
    for ins in data:
        ins_label['images'].append({
            'file_name': ins['image_name'],
            'height': ins['height'],
            'width': ins['width'],
            'id': image_id,
            'image_name': '_'.join(ins['image_name'].split('.')[0].split('_')[-3:])
        })
        assert len(ins['bbox']) == len(ins['segmentation'])
        for bbox, seg in zip(ins['bbox'], ins['segmentation']):
            ins_label['annotations'].append({
                'area': bbox[2]*bbox[3],
                'iscrowd': 0,
                'bbox': bbox,
                'image_id': image_id,
                'category_id': 1,
                'id': anno_id,
                'ignore': 0,
                'segmentation': seg
            })
            anno_id += 1
        image_id += 1
    with open('./json/nocrop_ins_label.json', 'w') as f:
        j.dump(ins_label, f)


def bbox2coco(data):
    ins_label = {
        'images': [],
        'type': 'instance',
        'annotations': [],
        'categories': [{'supercategory': 'none', 'id': 1, 'name': 'granule'}]
    }
    image_id = 1
    anno_id = 1
    for ins in data:
        ins_label['images'].append({
            'file_name': ins['image_name'],
            'height': ins['height'],
            'width': ins['width'],
            'id': image_id,
            'image_name': '_'.join(ins['image_name'].split('.')[0].split('_')[-3:])
        })
        assert len(ins['bbox']) == len(ins['segmentation'])
        for bbox, seg in zip(ins['bbox'], ins['segmentation']):
            ins_label['annotations'].append({
                'area': bbox[2]*bbox[3],
                'iscrowd': 0,
                'bbox': bbox,
                'image_id': image_id,
                'category_id': 1,
                'id': anno_id,
                'ignore': 0,
                'segmentation': seg
            })
            anno_id += 1
        image_id += 1
    with open('/p300/ihuman/debug/temp2.json', 'w') as f:
        j.dump(ins_label, f)
    return ins_label


if __name__ == "__main__":
    process()


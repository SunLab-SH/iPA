import os
import sys
import cv2 as cv
import json as j
from math import ceil
from time import sleep
import argparse as arg
import pycocotools.mask as mask_util
from itertools import groupby
from schedule import Schedule


def error(msg='Unknow'):
    print('\nERROR: ' + msg)
    exit(1)


def _fix_path(p):
    if p[-1] != '/':
        return p + '/'
    else:
        return p


def show_cfg(cfg):
    if cfg.search == 1:
        cfg.step = round(cfg.step, 2)
        print(f'>> Mode: search')
        if cfg.step <= 0 or cfg.step >= 1:
            error('step must be set to (0,1)!')
        print(f'>> Step: {cfg.step}')
    else:
        print(f'>> Mode: transform')

    if cfg.mode == 0:
        print('>> Single Image: only split one image')
        print(f'>> Image path: "{cfg.path}"')
        if cfg.ignore == 0 and not os.path.isfile(cfg.path):
            error('image not found!')
    else:
        print('Multiple Images: split images in directory')
        print(f'>> Image dir: "{cfg.dir}"')
        if cfg.ignore == 0 and not os.path.isdir(cfg.dir):
            error('dir not found!')
        cfg.dir = _fix_path(cfg.dir)

    if cfg.overlap == -1:
        print(f'>> Split Width: {cfg.width}')
        print(f'>> Split Height: {cfg.height}')
    else:
        print(f'>> Overlap: {cfg.overlap}')
    print(f'>> Split Grid: {cfg.xgrid}X{cfg.ygrid}')
    print(f'>> Top limit: {cfg.limit}')
    print(f'>> Annotaitons: "{cfg.anno}"')
    if cfg.ignore == 0 and not os.path.isfile(cfg.anno):
        error('annotation file not found!')
    print(f'>> Store dir: "{cfg.store}"')
    if cfg.ignore == 0 and not os.path.isdir(cfg.store):
        error('store dir not found!')
    cfg.store = _fix_path(cfg.store)
    if len(os.listdir(cfg.store)) != 0:
        if os.name == 'nt':
            new_dir = cfg.store.replace('/', '\\')
            # os.system(('del %strain\*.png')%(new_dir))
        else:
            pass
            # os.system(('rm %s')%(os.path.join(cfg.store, 'train/*.png')))
    if cfg.imgcrt == 0 and cfg.annocrt == 0:
        error(f'You must set a criterion of images list generation')
    if cfg.imgcrt != 0 and cfg.annocrt != 0:
        error(f'You can not set imgcrt and annocrt together')
    if cfg.imgcrt != 0:
        print(f'Use imgcrt')
    if cfg.annocrt != 0:
        print(f'Use annocrt')
    print('----------------------')


def store(cfg, new_anno_dict, result, img_id, anno_id, global_split_info):
    outdir = cfg.store
    for bimg in result:
        for origin_id, new_id, simg, new_bbox, new_seg, split_info in bimg:
            if len(new_bbox) < cfg.min_bbox:
                continue
            img_name = f'split_{origin_id}_{img_id}.png'
            global_split_info[img_name] = split_info
            # TODO: 修改outdir
            # name = '%s/train/%s' % (outdir, img_name)
            name = '/p300/ihuman/dataset/split_images/train/%s' % img_name
            # cv.imwrite(name, simg)
            h, w, _ = simg.shape
            new_anno_dict['images'].append({
                'file_name': img_name,
                'height': h,
                'width': w,
                'id': img_id,
                'image_name': '%d_%d' % (origin_id, img_id)
            })
            for bbox, seg in zip(new_bbox, new_seg):
                new_anno_dict['annotations'].append({
                    'area': bbox[2] * bbox[3],
                    'iscrowd': 0,
                    'image_id': img_id,
                    'bbox': bbox,
                    'category_id': 1,
                    'segmentation': seg,
                    'id': anno_id
                })
                anno_id += 1
            img_id += 1
    return new_anno_dict, img_id, anno_id


def get_anno(path):
    print('>> Loading Annotations...')
    f = open(path, 'r')
    d = j.load(f)
    print('<< Successfully Loaded!')
    print('>> Optimization processing...')
    img2id = {}
    for img in d['images']:
        img2id[img['file_name'].lower()] = img['id']
    id2bbox = {}
    id2seg = {}
    for anno in d['annotations']:
        if anno['image_id'] not in id2bbox:
            id2bbox[anno['image_id']] = [anno['bbox']]
            id2seg[anno['image_id']] = [anno['segmentation']]
        else:
            id2bbox[anno['image_id']].append(anno['bbox'])
            id2seg[anno['image_id']].append(anno['segmentation'])
    print('<< Successfully optimized!')

    # check number of bbox and seg is equal
    try:
        error_id = None
        for id in id2bbox.keys():
            error_id = id
            assert len(id2bbox[id]) == len(id2seg[id])
    except AssertionError as e:
        print(f'bbox != seg  id={error_id}')
        exit(1)
    print(f'Anno check successfully!')

    return (img2id, id2bbox, id2seg)


def get_original_id(anno, path):
    if '/' in path:
        path = path.split('/')[-1]
    name = path
    for img in anno['images']:
        if img['file_name'] == name:
            return img['id']
    error(('original id not found! img-path: %s') % (path))


def get_original_bbox(anno, id):
    results = []
    for box in anno['annotations']:
        if id == box['image_id']:
            results.append(box['bbox'])
    return results


def _calc_wh(raw, grid, overlap):
    return (2 * raw) / (2 * grid + overlap - overlap * grid)


def fix_avgwh(cfg, raw):
    best_grid = -1
    avg = cfg.avg
    min_gap = sys.maxsize
    for grid in range(1, 21):
        now_res = _calc_wh(raw, grid, cfg.overlap)
        if abs(now_res - avg) < min_gap:
            min_gap = abs(now_res - avg)
            best_grid = grid
    return best_grid


def fix_wh(cfg, raw_w, raw_h):
    if cfg.overlap == -1:
        w = int(max(cfg.width, ceil(raw_w / cfg.xgrid)))
        h = int(max(cfg.height, ceil(raw_h / cfg.ygrid)))
        # if w != cfg.width:
        #     print(f'>> No overlap: Width too small and fix width to {w}')
        # if h != cfg.height:
        #     print(f'>> No overlap: height too small and fix height to {h}')
    else:
        w = int((2 * raw_w) / (2 * cfg.xgrid + cfg.overlap - cfg.overlap * cfg.xgrid))
        h = int((2 * raw_h) / (2 * cfg.ygrid + cfg.overlap - cfg.overlap * cfg.ygrid))
        # if w != cfg.width:
        #     print(f'>> Overlap not satisfy: fix width to {w}')
        # if h != cfg.height:
        #     print(f'>> Overlap not satisfy: fix height to {h}')
    return (w, h)


def _fix_point(minv, maxv, origin, avoid):
    if origin < minv:
        if avoid == -1:
            error('bbox not in subimg')
        else:
            return minv
    elif origin > maxv:
        if avoid == 1:
            error('bbox not in subimg')
        else:
            return maxv
    else:
        return origin


def conflict_retain_large(threshold, x, y, w, h, bbox):
    bx, by, bw, bh = bbox

    new_x = _fix_point(x, x + w, bx, 1)
    new_y = _fix_point(y, y + h, by, 1)
    new_w = _fix_point(x, x + w, bx + bw, -1) - new_x
    new_h = _fix_point(y, y + h, by + bh, -1) - new_y
    new_x = new_x - x
    new_y = new_y - y

    raw_area = bw * bh
    new_area = new_w * new_h

    if new_area / raw_area < 0.5:
        return (False, ())
    else:
        return (True, (new_x, new_y, new_w, new_h))


def new_retain_large(threshold, raw_bbox, new_bbox):
    raw_area = raw_bbox[2] * raw_bbox[3]
    new_area = new_bbox[2] * new_bbox[3]
    if new_area / raw_area < threshold or new_area < 4:
        return False
    else:
        return True


def compute_bbox(cfg, x, y, w, h, bboxs):
    result = []
    conflict = 0
    w -= 1
    h -= 1
    seg_info = []
    for bbox in bboxs:
        bx, by, bw, bh = bbox
        now_conf = False
        new_bbox = []
        # 3 6 9 8 7
        if bx >= x + w or by >= y + h:
            seg_info.append(0)
            continue
        else:
            # 5
            if bx >= x and by >= y:
                new_x, new_y = bx - x, by - y
                if bx + bw > x + w:
                    new_w = x + w - bx
                    now_conf = True
                else:
                    new_w = bw
                if by + bh > y + h:
                    new_h = y + h - by
                    now_conf = True
                else:
                    new_h = bh
                new_bbox = [new_x, new_y, new_w, new_h]
            # 2
            if bx >= x and by < y:
                if by + bh <= y:
                    seg_info.append(0)
                    continue
                else:
                    now_conf = True
                    new_x, new_y = bx - x, 0
                    if bx + bw <= x + w:
                        new_w = bw
                    else:
                        new_w = x + w - bx
                    if by + bh <= y + h:
                        new_h = by + bh - y
                    else:
                        new_h = y
                    new_bbox = [new_x, new_y, new_w, new_h]
            # 4
            if bx < x and by >= y:
                if bx + bw <= x:
                    seg_info.append(0)
                    continue
                else:
                    now_conf = True
                    new_x, new_y = 0, by - y
                    if bx + bw <= x + w:
                        new_w = bx + bw - x
                    else:
                        new_w = w
                    if by + bh <= y + h:
                        new_h = bh
                    else:
                        new_h = y + h - by
                    new_bbox = [new_x, new_y, new_w, new_h]
            # 1
            if bx < x and by < y:
                if bx + bw < x or by + bh < y:
                    seg_info.append(0)
                    continue
                else:
                    now_conf = True
                    new_x, new_y = 0, 0
                    if bx + bw <= x + w:
                        new_w = bx + bw - x
                    else:
                        new_w = w
                    if by + bh <= y + h:
                        new_h = by + bh - y
                    else:
                        new_h = y
                    new_bbox = [new_x, new_y, new_w, new_h]

            if new_bbox == []:
                print(f'\nx:{x} y:{y} w:{w} h:{h}  bx:{bx} by:{by} bw:{bw} bh:{bh}')
                error('compute new bbox unknown situation')

            if now_conf:
                if new_retain_large(cfg.threshold, bbox, new_bbox):
                    seg_info.append(2)
                    result.append(new_bbox)
                    conflict += 1
                else:
                    seg_info.append(0)
            else:
                result.append(new_bbox)
                seg_info.append(1)

        """ old strategt (exist some bugs)
        # certainly not in subimg
        if bx > x+w or bx+bw < x or by > y+h or by+bh < y:
            continue
        else:
            # certainly in subimg
            if bx >= x and bx+bw <= x+w and by >= y and by+bh <= y+h:
                new_bbox = [bx-x, by-y, bw, bh]
                result.append(new_bbox)
            else:
                # certainly part of bbox in subimg
                if (bx >= x and bx <= x+w and by >= y and by <= y+h) or (bx >= x and bx <= x+w and by+bh >= y and by+bh <= y+h) or (bx+bw >= x and bx+bw <= x+w and by >= y and by <= y+h) or (bx+bw >= x and bx+bw <= x+w and by+bh >= y and by+bh <= y+h):
                    if cfg.strategy == 'pass':
                        conflict += 1
                    elif cfg.strategy == 'large':
                        retain,new_bbox = conflict_retain_large(cfg.threshold, x, y, w, h, bbox)
                        if retain:
                            result.append(new_bbox)
                            conflict += 1
                else:
                    print(f'\nx:{x} y:{y} w:{w} h:{h}  bx:{bx} by:{by} bw:{bw} bh:{bh}')
                    error('compute new bbox unknow situation')
        """
    assert len(seg_info) == len(bboxs)
    return (result, conflict, seg_info)


def binary_mask_to_rle(binary_mask):
    rle = {'counts': [], 'size': list(binary_mask.shape)}
    counts = rle.get('counts')
    for i, (value, elements) in enumerate(groupby(binary_mask.ravel(order='F'))):
        if i == 0 and value == 1:
            counts.append(0)
        counts.append(len(list(elements)))
    return rle


def compute_seg(x, y, w, h, segs, seg_info):
    result = []
    assert len(seg_info) == len(segs)
    for index, seg in enumerate(segs):
        isretain = seg_info[index]
        if isretain == 0:
            continue  # not retain bbox so there is no need to compute segmentation
        # print(f'seg: {seg}')
        compressed = mask_util.frPyObjects(seg, seg['size'][0], seg['size'][1])
        bmask = mask_util.decode(compressed)
        # although there is no difference but wait for later change, there retains if
        if isretain == 1:
            # there is not conflict between bbox and region
            bmask = bmask[y:y + h, x:x + w]
        elif isretain == 2:
            # there is conflict betwwen bbox and region
            bmask = bmask[y:y + h, x:x + w]

        new_seg_rle = binary_mask_to_rle(bmask)
        result.append(new_seg_rle)
    return result


def process_one_img(cfg, path, sid, img2id, id2bbox, id2seg):
    img = cv.imread(cfg.dir + path)
    raw_h, raw_w, _ = img.shape
    if '/' in path:
        path = path.split('/')[-1]
    path = path.lower()
    origin_id = img2id[path]
    origin_bbox = id2bbox[origin_id]
    origin_seg = id2seg[origin_id]

    if cfg.avg == -1:
        w, h = fix_wh(cfg, raw_w, raw_h)
    else:
        cfg.xgrid = fix_avgwh(cfg, raw_w)
        cfg.ygrid = fix_avgwh(cfg, raw_h)
        w, h = fix_wh(cfg, raw_w, raw_h)
    if cfg.xgrid != 1:
        x_step = (raw_w - w) / (cfg.xgrid - 1)
    else:
        x_step = 0
    x_points = [int(i * x_step) for i in range(cfg.xgrid)]
    if cfg.ygrid != 1:
        y_step = (raw_h - h) / (cfg.ygrid - 1)
    else:
        y_step = 0
    y_points = [int(i * y_step) for i in range(cfg.ygrid)]

    conflict = 0
    simgs = []
    for x in x_points:
        for y in y_points:
            split_info = (y, x)
            new_bbox, add_conflict, seg_info = compute_bbox(cfg, x, y, w, h, origin_bbox)
            conflict += add_conflict
            new_seg = compute_seg(x, y, w, h, origin_seg, seg_info)
            # new_seg = [1,1,2,1,2,2,1,1]
            assert len(new_bbox) == len(new_seg)

            simgs.append((origin_id, sid, img[y:y + h, x:x + w], new_bbox, new_seg, split_info))
            sid += 1

    return simgs, sid, conflict


def getpath_annocrt(cfg, img2id):
    real = os.listdir(cfg.dir)
    paths = real.copy()
    for i in range(len(paths)):
        paths[i] = paths[i].lower()
    result = []
    for key in img2id.keys():
        s = key.lower()
        if s in paths:
            pos = paths.index(s)
            result.append(real[pos])
        else:
            print(f'key: {key}')
            raise TypeError
    return result


def process(cfg, img2id, id2bbox, id2seg):
    sid = 1
    result = []
    new_anno_dict = {
        'annotations': [],
        'images': [],
        'type': 'instance',
        'categories': [{'supercategory': 'none', 'id': 1, 'name': 'granule'}]
    }
    BUFFER = 5
    conflict = 0
    global_split_info = {}

    if cfg.mode == 0:
        path = cfg.path
        res, sid, add_conflict = process_one_img(cfg, path, sid, img2id, id2bbox)
        conflict += add_conflict
        result.append(res)
        if cfg.search == 0:
            new_anno_dict = store(cfg, new_anno_dict, result, 0, 0)
    else:
        if cfg.imgcrt == 1:
            paths = os.listdir(cfg.dir)
        elif cfg.annocrt == 1:
            paths = getpath_annocrt(cfg, img2id)
        img_id = 1
        anno_id = 1
        if cfg.limit == -1:
            limit = len(paths)
        else:
            limit = cfg.limit
        if cfg.search == 0:
            print('')
        schedule = Schedule(limit)
        for idx in range(limit):
            path = paths[idx]
            res, sid, add_conflict = process_one_img(cfg, path, sid, img2id, id2bbox, id2seg)
            conflict += add_conflict
            result.append(res)
            # bar = '#' * (int(((idx + 1) / limit) * 40)) + '.' * (40 - int(((idx + 1) / limit) * 40))
            if cfg.search == 0:
                if len(result) >= BUFFER or idx == limit - 1:
                    new_anno_dict, img_id, anno_id = store(cfg, new_anno_dict, result, img_id, anno_id, global_split_info)
                    result = []
                # print('\rProcess: %d/%d [%s] %.2f%%' % (idx + 1, limit, bar, ((idx + 1) / limit) * 100), end='')
                schedule.watch()
            else:
                if idx > 100:
                    print('\rSubProcess: %d/%d [%s] %.2f%%' % (idx + 1, limit, bar, ((idx + 1) / limit) * 100),
                          end='')
    if cfg.search == 0:
        print('\n')

    # Store Split information
    with open(os.path.join(cfg.store, 'train', 'split_info_train.json'), 'w', encoding='utf-8') as w:
        j.dump(global_split_info, w)

    if cfg.search == 0:
        print('>> Storing annotations...')
        anno_out_name = os.path.join(cfg.store, 'train', 'split_ins_label_rle_train.json')
        with open(anno_out_name, 'w', encoding='utf-8') as f:
            j.dump(new_anno_dict, f)
        print(f'<< Successfully Stored! -> {anno_out_name}')
        print(
            f'>> Total images: {len(new_anno_dict["images"])}  Total annotations: {len(new_anno_dict["annotations"])}')
        print(f'>> Total conflict: {conflict}')
        return new_anno_dict
    else:
        return conflict


def compute_steps(step):
    result = []
    start = 0.0
    while start < 1.0:
        result.append(round(start, 2))
        start += step
    result.remove(0)
    return result


def compute_table_space(max_size, now):
    left = int((max_size - len(now)) / 2)
    right = max_size - left - len(now)
    # left += 1
    # right += 1
    return ' ' * left + now + ' ' * right


def show_search_result(perfect, min_conflict, all_conflict, steps):
    print('\n\n-------------------------')
    print('>> Search Result')
    for i in perfect:
        print(f'Perfect Overlap: {i}')
    print(f'Total conflict: {min_conflict}\n')

    max_size = -1
    for c in all_conflict:
        max_size = max(len(str(c)), max_size)
    for s in steps:
        max_size = max(len(str(s)), max_size)

    print('%s' % ('=' * int((max_size + 1) * len(steps) + 11)))
    print('| Overlap |', end='')
    line = '|---------|'
    for step in steps:
        display = compute_table_space(max_size, str(step))
        print(f'{display}|', end='')
        line += ('-' * max_size + '|')
    print('\n%s' % line)
    print('|Conflicts|', end='')
    for conflict in all_conflict:
        display = compute_table_space(max_size, str(conflict))
        print(f'{display}|', end='')
    print('\n%s' % ('=' * int((max_size + 1) * len(steps) + 11)))


"""
RLE -> Polygon 模块开始
"""

import sys
import numpy as np
import json as j
import pycocotools.mask as mask_util
from imantics import Polygons, Mask
from scipy.ndimage import binary_fill_holes, binary_dilation
from schedule import Schedule

"""
将标注文件中segmentation的RLE格式转换为Polygon格式，因为模型输入需要Polygon格式
因为转换方法非官方（官方未提供）因此转换后与RLE格式的结果会有差别
此步骤主要用于在split出小区域图片后对小区域图片标注文件进行转换，若用于别的转换也OK
"""


def rle2poly(data):
    schedule = Schedule(len(data['annotations']))
    for index in range(len(data['annotations'])):
        anno = data['annotations'][index]

        # segmentation
        rle = anno['segmentation']
        compressed_rle = mask_util.frPyObjects(rle, rle.get('size')[0], rle.get('size')[1])
        bmask = mask_util.decode(compressed_rle)
        bmask = binary_dilation(bmask).astype(int)
        polygons = Mask(bmask).polygons()
        data['annotations'][index]['segmentation'] = polygons.segmentation
        data['annotations'][index]['iscrowd'] = 0

        # bbox recalculate
        height = bmask.shape[0]
        width = bmask.shape[1]
        compressed = mask_util.frPyObjects(polygons.segmentation, height, width)
        new_bmask = mask_util.decode(compressed)
        points = np.where(bmask == 1)
        xs = points[0]
        ys = points[1]
        x_min = sys.maxsize
        y_min = sys.maxsize
        x_max = -1
        y_max = -1
        for x, y in zip(xs, ys):
            x_min = min(x, x_min)
            y_min = min(y, y_min)
            x_max = max(x, x_max)
            y_max = max(y, y_max)
        x = int(x_min) - 1
        y = int(y_min) - 1
        w = int(x_max + 1 - x_min + 1)
        h = int(y_max + 1 - y_min + 1)
        data['annotations'][index]['bbox'] = [y, x, h, w]
        schedule.watch()
    return data

"""
RLE -> Poly 模块结束
"""


def init():
    parser = arg.ArgumentParser(description='Split granule img into 20 small imgs')
    parser.add_argument('-m', '--mode', dest='mode', default=1, type=int)
    parser.add_argument('-p', '--path', dest='path', default='/root/mask-rcnn/split/special.png', type=str)
    parser.add_argument('-d', '--dir', dest='dir', default='/group/xiangyi/Chuanyang/denoise_images/train/', type=str)
    parser.add_argument('-l', '--limit', dest='limit', default=-1, type=int)
    parser.add_argument('-s', '--store', dest='store', default='/root/mask-rcnn/split_smaller/', type=str)
    parser.add_argument('-a', '--annotations', dest='anno',
                        default='/group/xiangyi/Chuanyang/annotations/crop_change_instances_train2014.json', type=str)
    parser.add_argument('--min-bbox', dest='min_bbox', default=3, help='留下的标注中最少的bbox数量')

    parser.add_argument('--width', dest='width', default=100, type=float)
    parser.add_argument('--height', dest='height', default=150, type=float)
    parser.add_argument('-o', '--overlap', dest='overlap', default=-1, type=float)
    parser.add_argument('-avg', '--average', dest='avg', default=-1, type=float)
    parser.add_argument('-x', '--xgrid', dest='xgrid', default=4, type=int)
    parser.add_argument('-y', '--ygrid', dest='ygrid', default=5, type=int)
    parser.add_argument('--ignore', dest='ignore', default=0, type=int)
    parser.add_argument('--search', dest='search', default=0, type=int)
    parser.add_argument('--step', dest='step', default=0.1, type=float)
    parser.add_argument('--strategy', dest='strategy', default='pass', type=str)
    parser.add_argument('--threshold', dest='threshold', default=0.5, type=float)
    parser.add_argument('--imgcrt', default=0, type=int, help='list images which will be splited by os.listdir')
    parser.add_argument('--annocrt', default=0, type=int, help='list images which will be splited by annotation file')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    cfg = init()
    show_cfg(cfg)
    img2id, id2bbox, id2seg = get_anno(cfg.anno)

    if cfg.search == 0:
        rle_label = process(cfg, img2id, id2bbox, id2seg)
        poly_out_path = os.path.join(cfg.store, 'train', 'split_ins_label_poly_train.json')
        print('Generate polygon annotation file:')
        poly_label = rle2poly(rle_label)
        with open(poly_out_path, 'w') as f:
            j.dump(poly_label, f)
    else:
        perfect = []
        all_conflict = []
        min_conflict = sys.maxsize
        steps = compute_steps(cfg.step)
        print('')
        for index, overlap in enumerate(steps):
            if index == 0:
                bar = '#' * (int((index + 1) / len(steps) * 20)) + '.' * (20 - int((index + 1) / len(steps) * 20))
                print('\rProcess: %d/%d [%s] %.2f%%' % (index + 1, len(steps), bar, (index + 1) / len(steps) * 100),
                      end='')
            cfg.overlap = overlap
            now_conflict = process(cfg, img2id, id2bbox, id2seg)
            all_conflict.append(now_conflict)
            if now_conflict == min_conflict:
                perfect.append(overlap)
            elif now_conflict < min_conflict:
                perfect = [overlap]
                min_conflict = now_conflict
            bar = '#' * (int((index + 1) / len(steps) * 20)) + '.' * (20 - int((index + 1) / len(steps) * 20))
            print('\rMainProcess: %d/%d [%s] %.2f%% %s' % (
                index + 1, len(steps), bar, (index + 1) / len(steps) * 100, ' ' * 30), end='')
        show_search_result(perfect, min_conflict, all_conflict, steps)





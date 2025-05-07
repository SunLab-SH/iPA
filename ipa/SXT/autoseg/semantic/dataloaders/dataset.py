# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 23:09
# @Author  : Xiangyi Zhang
# @File    : dataset.py
# @Email   : zhangxy9@shanghaitech.edu.cn

from cgi import test
import cv2
import numpy as np
import os
import sys
import glob
import itertools
import torch
from torch.utils.data import Dataset
import dataloaders.transforms as tf
from torch.utils.data import DataLoader
from torchvision.transforms import Compose, ColorJitter
from torch.utils.data.sampler import Sampler
from experiments.parser import args

sys.path.append(os.path.abspath('./'))
image_size = (288, 480)
# width, height


class ImageLoader(Dataset):
    def __init__(self, mode, data_root_dir=None, transforms=None):
        self.mode = mode
        self.all_label = True
        self.channel3 = True
        self.transforms = transforms
        assert data_root_dir is not None, "Please pass the correct the data_root_dir"
        assert self.mode in ["train", "val", "test", "test_refine"], "Only support train, val, test and unlabel"

        image_dir = 'images_file/'
        label_dir = 'labels_file/'
        images_list = []
        for path in os.listdir(os.path.join(data_root_dir, image_dir, self.mode)):
            images_list.extend(glob.glob(os.path.join(data_root_dir, image_dir, self.mode, path)))
        images_list.sort()

        labels_list = []
        for path in os.listdir(os.path.join(data_root_dir, image_dir, self.mode)):
            labels_list.extend(glob.glob(os.path.join(data_root_dir, label_dir, self.mode, path)))
        labels_list.sort()

        # print(len(labels_list), len(images_list))
        assert (len(labels_list) == len(images_list)), "lens of labels and images should be same."

        self.images_list = images_list
        self.labels_list = labels_list

        print('Done initializing ' + self.mode + ' Dataset')

    def __getitem__(self, item):
        image = cv2.imread(os.path.join(self.images_list[item]), 0)  # X, Y, 3
        image = cv2.resize(image, image_size, cv2.INTER_NEAREST)
        label = cv2.imread(os.path.join(self.labels_list[item]), 0)
        label = cv2.resize(label, image_size, cv2.INTER_NEAREST)

        if args.num_classes == 2:
            label[label == 1] = 2
            label[label == 4] = 2
            label[label == 2] = 1
            label[label == 3] = 2
        # There is no need to do transpose before the transform,
        # it will do transform in tf.Totensor()
        sample = {'image': image, 'label': label}
        if self.transforms is not None:
            sample = self.transforms(sample)
        else:
            sample["image"] = np.transpose(sample["image"], axes=[2, 0, 1])
        return sample

    def __len__(self):
        return len(self.images_list)


class UnlabelImageDataloader(Dataset):
    def __init__(self, data_root_dir=None, transforms=None):
        """
        load unlabel images which do not have segmentation label
        :param data_root_dir: the root path to the dataset, and default unlabel image dir is "unlabel"
        :return: None
        """
        self.transforms = transforms
        images_list = []
        if args.test_idx == 'all':
            for path in os.listdir(os.path.join(data_root_dir)):
                images_list.extend(glob.glob(os.path.join(data_root_dir, path, "*.png")))
        elif args.test_idx == 'iso':
            images_list.extend(glob.glob(os.path.join(data_root_dir, "*.png")))
        else:
            images_list.extend(glob.glob(os.path.join(data_root_dir, args.test_idx, "*.png")))
        images_list.sort()
        self.images_list = images_list

    def __getitem__(self, item):

        image = cv2.imread(os.path.join(self.images_list[item]), 0)  # X, Y, 3
        image = cv2.resize(image, image_size, cv2.INTER_NEAREST)

        sample = {'image': image}
        if self.transforms is not None:
            sample = self.transforms(sample)
        return sample

    def __len__(self):
        return len(self.images_list)


class ImageLoaderLabelUnlabel(Dataset):
    def __init__(self, data_root_dir=None, transforms=None):
        self.all_label = True
        self.channel3 = True
        self.transforms = transforms
        assert data_root_dir is not None, "Please pass the correct the data_root_dir"

        image_dir = 'images_file/'
        label_dir = 'labels_file/'
        images_list = []
        for path in os.listdir(os.path.join(data_root_dir, image_dir, "train")):
            images_list.extend(glob.glob(os.path.join(data_root_dir, image_dir, "train", path, "*.png")))
        images_list.sort()

        labels_list = []
        for path in os.listdir(os.path.join(data_root_dir, image_dir, "train")):
            labels_list.extend(glob.glob(os.path.join(data_root_dir, label_dir, "train", path, "*.png")))
        labels_list.sort()
        assert len(labels_list) == len(images_list), "Num of label files and num of image files is different."

        unlabel_images_list = []
        for path in os.listdir(os.path.join(data_root_dir, image_dir, "unlabel")):
            unlabel_images_list.extend(glob.glob(os.path.join(data_root_dir, image_dir, "unlabel", path, "*.png")))
        unlabel_images_list.sort()

        # print(len(labels_list), len(images_list))
        self.label_lens = len(images_list)
        self.unlabel_lens = len(unlabel_images_list)

        images_list.extend(unlabel_images_list)
        self.images_list = images_list
        self.labels_list = labels_list

        print('Done initializing ' + "label and unlabel" + ' Dataset')

    def __getitem__(self, item):
        image = cv2.imread(os.path.join(self.images_list[item]), 0)  # X, Y, 3
        # Image size after resize is (480, 288)
        image = cv2.resize(image, image_size, cv2.INTER_NEAREST)
        try:
            label = cv2.imread(os.path.join(self.labels_list[item]), 0)
            label = cv2.resize(label, image_size, cv2.INTER_NEAREST)
        except IndexError:
            label = np.zeros((image_size[1], image_size[0]))

        if not self.all_label:
            label[label == 2] = 0
            label[label == 3] = 0
            label[label == 4] = 2
        # There is no need to do transpose before the transform,
        # it will do transform in tf.Totensor()
        sample = {'image': image, 'label': label}
        if self.transforms is not None:
            sample = self.transforms(sample)
        else:
            sample["image"] = np.transpose(sample["image"], axes=[2, 0, 1])
        return sample

    def __len__(self):
        return len(self.images_list)


def data_loader(paras):
    train_transform = Compose([tf.RandomHorizontalFlip(),
                               tf.ScaleNRotate(rots=(-10, 10), scales=(.8, 1.2)),
                               # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                               tf.ToTensor(),
                               tf.Normalize(mean=[0.456], std=[0.224])])

    val_transform = Compose([tf.ToTensor(),
                             tf.Normalize(mean=[0.456], std=[0.224])])

    train_dataset = ImageLoader(mode='train', data_root_dir=paras.data_root_dir,
                                transforms=train_transform)

    train_loader = DataLoader(train_dataset, batch_size=paras.batch_size,
                              shuffle=True, num_workers=paras.num_workers)

    val_dataset = ImageLoader(mode='val', data_root_dir=paras.data_root_dir,
                              transforms=val_transform)

    val_loader = DataLoader(val_dataset, batch_size=paras.batch_size,
                            shuffle=False, num_workers=paras.num_workers)

    test_dataset = ImageLoader(mode='test', data_root_dir=paras.data_root_dir,
                               transforms=val_transform)

    test_loader = DataLoader(test_dataset, batch_size=paras.batch_size,
                             shuffle=False, num_workers=paras.num_workers)

    return train_loader, val_loader, test_loader


def get_data_fs(args):
    """
    Dataloader for fully supervised experiment.
    """
    train_transform = Compose([tf.RandomHorizontalFlip(),
                               tf.ScaleNRotate(rots=(-10, 10), scales=(.8, 1.2)),
                               # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                               tf.ToTensor(),
                               tf.Normalize(mean=[0.456], std=[0.224])])

    val_transform = Compose([tf.ToTensor(),
                             tf.Normalize(mean=[0.456], std=[0.224])])

    train_dataset = ImageLoader(mode='train', data_root_dir=args.data_root_dir,
                                transforms=train_transform)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size,
                              shuffle=True, num_workers=args.num_workers)

    val_dataset = ImageLoader(mode='val', data_root_dir=args.data_root_dir,
                              transforms=val_transform)

    val_loader = DataLoader(val_dataset, batch_size=args.batch_size,
                            shuffle=False, num_workers=args.num_workers)

    return train_loader, val_loader


def get_data_test(args):
    val_transform = Compose([tf.ToTensor(),
                             tf.Normalize(mean=[0.456], std=[0.224])])

    val_dataset = ImageLoader(mode='val', data_root_dir=args.data_root_dir,
                              transforms=val_transform)

    val_loader = DataLoader(val_dataset, batch_size=args.batch_size,
                            shuffle=False, num_workers=args.num_workers)

    test_dataset = ImageLoader(mode='test', data_root_dir=args.data_root_dir,
                               transforms=val_transform)

    test_loader = DataLoader(test_dataset, batch_size=args.batch_size,
                             shuffle=False, num_workers=args.num_workers)
    return val_loader, test_loader


def get_data_un(args, data_root_dir):
    """
    Dataloader for semi supervised experiment.
    """
    transform = Compose([tf.ToTensor(),
                         tf.Normalize(mean=[0.456], std=[0.224]),
                         # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                         ])
    unlabel_dataset = UnlabelImageDataloader(data_root_dir=data_root_dir, transforms=transform)
    # test_imgs = unlabel_dataset['image']
    # # temp_img = images[i].cpu().numpy()
    # print(test_imgs.shape)
    
    unlabel_loader = DataLoader(unlabel_dataset, batch_size=args.batch_size, shuffle=False,
                                      num_workers=args.num_workers)
    # unlabel_loader = unlabel_dataset
    return  unlabel_loader


def get_data_un2(args, data_root_dir):
    """
    Dataloader for semi supervised experiment.
    """
    transform = Compose([tf.ToTensor(),
                         tf.Normalize(mean=[0.456], std=[0.224]),
                         # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                         ])
    unlabel_dataset = UnlabelImageDataloader(data_root_dir=data_root_dir, transforms=transform)
    # test_imgs = unlabel_dataset
    # temp_img = images[i].cpu().numpy()
    # print(test_imgs.shape)
    # print(test_imgs)
    
    
    unlabel_loader = DataLoader(unlabel_dataset, batch_size=args.batch_size, shuffle=False,
                                      num_workers=args.num_workers)
    # unlabel_loader = unlabel_dataset
    return  unlabel_loader


def get_data_ss(args):
    """
    Dataloader for semi supervised experiment.
    """
    train_transform = Compose([tf.RandomHorizontalFlip(),
                               tf.ScaleNRotate(rots=(-10, 10), scales=(.8, 1.2)),
                               # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                               tf.ToTensor(),
                               tf.Normalize(mean=[0.456], std=[0.224])])

    val_transform = Compose([tf.ToTensor(),
                             tf.Normalize(mean=[0.456], std=[0.224])])

    train_dataset = ImageLoader(mode='train', data_root_dir=args.data_root_dir,
                                transforms=train_transform)

    train_loader = DataLoader(train_dataset, batch_size=args.label_batch,
                              shuffle=True, num_workers=args.num_workers)

    train_unlabel_dataset = UnlabelImageDataloader(data_root_dir=args.data_root_dir, transforms=train_transform)

    train_unlabel_loader = DataLoader(train_unlabel_dataset, batch_size=args.batch_size - args.label_batch, shuffle=True,
                                      num_workers=args.num_workers)

    val_dataset = ImageLoader(mode='val', data_root_dir=args.data_root_dir,
                              transforms=val_transform)

    val_loader = DataLoader(val_dataset, batch_size=args.batch_size,
                            shuffle=False, num_workers=args.num_workers)

    return train_loader, train_unlabel_loader, val_loader


def get_data_ss_dev(args):
    """
    Dataloader for semi supervised experiment.
    """
    train_transform = Compose([tf.RandomHorizontalFlip(),
                               tf.ScaleNRotate(rots=(-10, 10), scales=(.8, 1.2)),
                               # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
                               tf.ToTensor(),
                               tf.Normalize(mean=[0.456], std=[0.224])])

    val_transform = Compose([tf.ToTensor(),
                             tf.Normalize(mean=[0.456], std=[0.224])])

    # RuntimeError: invalid argument 0: sizes of tensor must match except in dimension 0. Got 480 288 in dimension 2.
    # Bug Fixed. Change np.zeros(image_size) to np.zeros((image_size[1], image_size[0])) in Line 151.
    # cv2 image size is (H, W, C), resize function parameters is (W, H). This is a bad design.
    # This function apply a different way combine label and unlabel dataloader. It will be slow because of try..except..
    train_dataset = ImageLoaderLabelUnlabel(data_root_dir=args.data_root_dir,
                                            transforms=train_transform)
    label_lens = train_dataset.label_lens
    unlabel_lens = train_dataset.unlabel_lens
    labeled_idxs = list(range(label_lens))
    unlabeled_idxs = list(range(label_lens, label_lens + unlabel_lens))
    batch_sampler = TwoStreamBatchSampler(labeled_idxs, unlabeled_idxs, args.batch_size, args.batch_size - args.label_batch)

    train_loader = DataLoader(train_dataset, batch_sampler=batch_sampler, num_workers=args.num_workers)

    val_dataset = ImageLoader(mode='val', data_root_dir=args.data_root_dir,
                              transforms=val_transform)

    val_loader = DataLoader(val_dataset, batch_size=args.batch_size,
                            shuffle=False, num_workers=args.num_workers)

    return train_loader, val_loader


class TwoStreamBatchSampler(Sampler):
    """Iterate two sets of indices
    An 'epoch' is one iteration through the primary indices.
    During the epoch, the secondary indices are iterated through
    as many times as needed.
    """
    def __init__(self, primary_indices, secondary_indices, batch_size, secondary_batch_size):
        self.primary_indices = primary_indices
        self.secondary_indices = secondary_indices
        self.secondary_batch_size = secondary_batch_size
        self.primary_batch_size = batch_size - secondary_batch_size

        assert len(self.primary_indices) >= self.primary_batch_size > 0
        assert len(self.secondary_indices) >= self.secondary_batch_size > 0

    def __iter__(self):
        primary_iter = iterate_once(self.primary_indices)
        secondary_iter = iterate_eternally(self.secondary_indices)
        return (
            primary_batch + secondary_batch
            for (primary_batch, secondary_batch)
            in zip(grouper(primary_iter, self.primary_batch_size), grouper(secondary_iter, self.secondary_batch_size))
        )

    def __len__(self):
        return len(self.primary_indices) // self.primary_batch_size


def iterate_once(iterable):
    return np.random.permutation(iterable)


def iterate_eternally(indices):
    def infinite_shuffles():
        while True:
            yield np.random.permutation(indices)
    return itertools.chain.from_iterable(infinite_shuffles())


def grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF"
    args = [iter(iterable)] * n
    return zip(*args)


def TwoDataLoader(label, unlabel):
    for lb_sample, ub_sample in zip(label, unlabel):
        lb_sample["image"] = torch.cat((lb_sample["image"], ub_sample["image"]), dim=0)
        yield lb_sample


if __name__ == "__main__":
    pass

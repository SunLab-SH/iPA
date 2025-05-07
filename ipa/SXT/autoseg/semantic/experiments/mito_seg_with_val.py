#%%
# -*- coding: utf-8 -*-
# @Time    : 2020-06-02 22:10
# @Author  : Xiangyi Zhang
# @File    : eval_mito.py
# @Email   : zhangxy9@shanghaitech.edu.cn



import os, sys
from unittest import TestLoader
filepath = os.path.dirname(os.path.abspath(__file__))
# print(filepath)
os.chdir(filepath)
print(os.getcwd())

os.chdir(os.path.pardir)
parpath_ = os.getcwd()
sys.path.append(parpath_)


from dataloaders.dataset import *
from experiments.parser_mito import args


from networks.Unet import Unet
from torch import nn
import logging
import glob
import cv2
LOG = logging.getLogger('main')
import os.path as osp
import tifffile, mrcfile
from utils.src import *

os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
from pylab import plt
from multiprocessing.dummy import Pool
import time


print(args)
# print(type(args.num_classes))
image_size = (288, 480) # (288, 480)

#%%
'''
1: mem
2: mito
3: nu
4: gr
test_idx = ['784_5', '766_8', '842_17']
val_idx = ['783_5', '766_5', '842_12']
train_idx = ['766_2', '766_7', '766_10', '766_11', '769_5', '769_7', '783_6', '783_12', '784_4', '784_6', '784_7', '785_7', '822_4', '822_6', '822_7', '842_13', '931_11', '931_14']
'''
def normalize(image):
    '''
    for 2d image to normalize 
    '''
    image = (image - np.min(image) )/ np.max(image) * 255
    return image.astype(np.uint8)


class load_data():
    def __init__(self, dir, dataid) -> None:
        self.dir = dir
        self.dataid = dataid

    def get_filelist(self, dir, Filelist):
        newDir = dir
        if os.path.isfile(dir):
            Filelist.append(dir)
            # Filelist.append(os.path.basename(dir))
        elif os.path.isdir(dir):
            for s in os.listdir(dir):
                # 如果需要忽略某些文件夹，使用以下代码
                #if s == "xxx":
                    #continue
                newDir=os.path.join(dir,s)
                self.get_filelist(newDir, Filelist)

        return Filelist


    def get_rawimgname(self):
        namelist = self.get_filelist(self.dir, [])
        # print(namelist)
        for name in namelist:
            if f'{self.dataid}' in name and 'pre_rec.tif' in name:
                return name
        else:
            print(' raw img not exist')
            os._exit(0)


    def get_labelname(self):
        namelist = self.get_filelist(self.dir,[])
        for name in namelist:
            if f'{self.dataid}' in name and 'membrane.tif' in name.lower():
                memmaskname = name
            elif f'{self.dataid}' in name and 'nucleus.tif' in name.lower():
                ncmaskname = name
            elif f'{self.dataid}' in name and 'mito.tif' in name.lower():
                mitomaskname = name
            elif f'{self.dataid}' in name and 'vesicles.tif' in name.lower():
                vesiclemaskname = name

        # else:
        #     print('mask not exist')
            # os._exit(0)
        # print([memmaskname, ncmaskname, mitomaskname, vesiclemaskname])
        return [memmaskname, ncmaskname, mitomaskname, vesiclemaskname]


    def get_labellfile(self):
        labelfilename = self.get_labelname()
        memimg = tifffile.imread(labelfilename[0])
        mitoimg = tifffile.imread(labelfilename[2])
        nucleusimg = tifffile.imread(labelfilename[1])
        vesicleimg = tifffile.imread(labelfilename[3])
        mask = np.zeros_like(memimg)
        mask[np.where(memimg != 0)] = 1
        mask[np.where(mitoimg != 0)] = 2
        mask[np.where(nucleusimg != 0)] = 3
        mask[np.where(vesicleimg != 0)] = 4

        return mask


# with open(os.path.join(args.data_root_dir, 'dataid.txt'), 'r') as f1:
#     data_idx = f1.read().splitlines()


# for id in data_idx:
#     data_a = load_data(os.path.join(args.data_root_dir), id)
#     print(data_a.get_rawimgname())
#     image_x = tifffile.imread(data_a.get_rawimgname())
#     print(image_x.shape)
#     print(image_x, np.max(image_x), np.min(image_x))
# #%%

# plt.imshow(image_x[234])
# plt.show()
# plt.close()

# y = np.rot90(image_x, k=3, axes=(2, 0))  ## y
# z = np.rot90(image_x, k=1, axes=(1, 0))  ## z
    
# plt.imshow(y[123])
# plt.show()
# plt.close()

# plt.imshow(z[234])
# plt.show()
# plt.close()


# #%%

# imagetest = image_x[234]
# imagetest = normalize(imagetest)

# image_size = (288, 480) # (288, 480)
# imagetest = cv2.resize(imagetest, image_size, cv2.INTER_NEAREST)
# plt.imshow(imagetest)
# plt.show()
# plt.close()

# sample = {'image': imagetest}

# transform = Compose([tf.ToTensor(),
#                         tf.Normalize(mean=[0.456], std=[0.224]),
#                         # tf.Normalize(mean=[0.456], std=[0.224]),
#                         # ColorJitter(brightness=0.8, contrast=0.8, saturation=0.4, hue=0.1),
#                         ])


# imagedata =transform(sample)
# print(np.max(imagedata['image'].cpu().numpy()))
# # print(imagedata['image'])


# image22 = imagedata['image']

# image22 = torch.unsqueeze(image22, 0)
# # image22 = torch.transpose(image22, -1, -2)
# print('line122', image22.shape)

# num_classes = 2
# checkpoint_dir = f'F:/salilab/salilab/projects/auto-segmentation_Aneesh/fromgit/Cell-Segmentation/results/FS_mito/best.pth'


# model = Unet(n_class=num_classes+1, is_dropout=True)
# # model = model.to(DEVICE)
# # model = nn.DataParallel(model)

# model = torch.nn.DataParallel(model)

# model.load_state_dict(torch.load(checkpoint_dir))
# model = model.module.to(torch.device('cpu'))

# # print(unlabel_loader)

# outputs = model(image22)
# pred = torch.argmax(torch.softmax(outputs, dim=1), dim=1)


# predd = torch.squeeze(pred, 0)
# predd = predd.cpu().numpy()
# plt.imshow(predd)
# # plt.title(f'{b}')
# plt.show()
# plt.close()

# #%%







#%%



class cellmapping_test():
    def __init__(self, type='imageall_x'):
        self.type = type
        self.image_root_dir = osp.join('data/image_xyz/', self.type)
        self.checkpoint_dir =  'F:/salilab/salilab/projects/auto-segmentation_Aneesh/fromgit/Cell-Segmentation/results/FS_mito/best.pth' # 'results/{}/best.pth'.format(args.exp)
        # self.device = torch.device("cuda")
        self.device = torch.device("cuda")
        self.gt_path = osp.join('data/mask_xyz')
        self.coarse_label_dir = 'data/coarse_label'

        # with open('test/test_idx.txt', 'r') as f1, open('test/val_idx.txt', 'r') as f2, open('test/train_idx.txt', 'r') as f3, open('test/unlabel_idx.txt', 'r') as f4:
        #     self.test_idx = f1.read().splitlines()
        #     self.val_idx = f2.read().splitlines()
        #     self.train_idx = f3.read().splitlines()
        #     self.unlabel_idx = f4.read().splitlines()
        with open(os.path.join(args.data_root_dir, 'dataid.txt'), 'r') as f1:
            data_idxs = f1.read().splitlines()
        self.data_idxs = data_idxs
        self.test_idx = self.data_idxs


    def save_tiff(self, pred, imagefile):
        # imagefile = imageid
        savepath = f'{args.root_dir}//results//masks//{self.type}_mito_predtiff'
        os.makedirs(savepath, exist_ok=True)
        
        tifffile.imsave(f'{savepath}/{imagefile}.tiff', np.array(pred))
        print(f"save {savepath}/{imagefile}.tiff'")

    def evaluate_iou(self, x, y, num_class):
        """
        :param x: tensor (B, H, W) {0, 5} float
        :param y: tensor (B, 1, H, W) {0, 5} float
        :return: IOU: float
        """
        batch_size, h, w = x.size()
        x = x.reshape(batch_size, -1).long()
        y = y.reshape(batch_size, -1).long()
        acc = accuracy(x, y)
        IOU, Dice = IOU_dice(x, y, num_class)
        return IOU, acc, Dice

    def resize(self, pred, image_ori_path=None):
        label_ori = cv2.imread(image_ori_path, 0)
        # plt.imshow(pred)
        # plt.show()
        # plt.close()
        pred_resize = cv2.resize(pred.astype(np.uint8), (label_ori.shape[1], label_ori.shape[0]), cv2.INTER_NEAREST)
        return pred_resize

    def get_path(self, image_root_dir):
        images_list = []

        if args.test_idx == 'all':
            for path in os.listdir(os.path.join(image_root_dir)):
                images_list.extend(glob.glob(os.path.join(image_root_dir, path, "*.png")))
        elif args.test_idx == 'iso':
            images_list.extend(glob.glob(os.path.join(image_root_dir, "*.png")))
        else:
            images_list.extend(glob.glob(os.path.join(image_root_dir, args.test_idx, "*.png")))
        images_list.sort()
        # print('line104', images_list)
        return images_list


    def turnGtBigLabel(self, imagefile):
        ##### Evaluate on membrane and nuclear
        ##### Set mito and nuclear as membrane(1), and nuclear(2)
        # gt = tifffile.imread(osp.join(self.gt_path, imagefile +'_merged_4organelles_mask.tiff'))
        data_a = load_data(os.path.join(args.data_root_dir), imagefile)
        gt = data_a.get_labellfile()
        # print('line110',gt.shape)
        if args.mode == 'mito':
        ## Small label: mito
            gt[gt != 2] = 0
            gt[gt == 2] = 1
        ## Big label: Membrane and Nuclear
        else:
            gt[gt == 1] = 2
            gt[gt == 4] = 2
            gt[gt == 2] = 1
            gt[gt == 3] = 2
        return gt

    def turnPredMitoLabel(self, pred):
        pred[pred != 1] = 0
        return pred

    def evaluate(self, pred_list, gt):
        pred_list = torch.Tensor(pred_list).cuda()
        gt = torch.Tensor(gt).cuda()
        # print('line129',pred_list)
        pred_list = pred_list[70:380, :,70:380]
        gt = gt[70:380, :,70:380]
        iou, acc, dice = self.evaluate_iou(pred_list, gt, num_class=args.num_classes + 1)
        return dice, iou, acc

    ### Inference ###
    def test(self):
        self.model = Unet(n_class=args.num_classes + 1, is_dropout=True)
        self.model = nn.DataParallel(self.model).cuda()
        self.model.load_state_dict(torch.load(self.checkpoint_dir))
        self.model = self.model.module.to(torch.device('cpu'))  # run on cpu
        print('Resume from {}'.format(self.checkpoint_dir))  # results/FS_mito/best.pth
        # print('line142',self.test_idx)


        print(self.data_idxs)
        for dataid in self.data_idxs:
            
            data_a = load_data(os.path.join(args.data_root_dir), dataid)
            print(data_a.get_rawimgname())
            image_x = tifffile.imread(data_a.get_rawimgname())

            # print(image_x.shape)
            # print(image_x, np.max(image_x), np.min(image_x))

            if self.type == 'imageall_x':
                image = image_x
            elif self.type == 'imageall_y':
                image = np.rot90(image_x, k=3, axes=(2, 0))
            elif self.type == 'imageall_z':
                image = np.rot90(image_x, k=1, axes=(1, 0))

            pred_mask = np.zeros_like(image)
            
            self.model.eval()

            def slicepred(ii):
                imageslice = image[ii]
                shape_ = imageslice.shape
                # print(shape_)
                
                inputslice = normalize(imageslice)
                inputslice = cv2.resize(inputslice, image_size, cv2.INTER_NEAREST)
                
                transform = Compose([tf.ToTensor(),
                        tf.Normalize(mean=[0.456], std=[0.224]),
                        ])
                sample = {'image': inputslice}
                imagedata =transform(sample)
                imageslice = imagedata['image']
                predslice = self.step(imageslice, dataid)
                predslice = predslice.cpu().numpy()
                # plt.imshow(predslice)
                # plt.show()
                # plt.close()
                predslice = cv2.resize(predslice.astype(np.uint8), (shape_[1], shape_[0]), cv2.INTER_NEAREST)
                
                pred_mask[ii] = predslice
                # time.sleep(1)

                if ii %10 == 0:
                    print(f'curii {ii}')
                    # plt.imshow(image[ii])
                    # plt.title(f'{ii}')
                    # plt.show()
                    # plt.close()
                    # plt.imshow(pred_mask[ii])
                    # plt.show()
                    # plt.close()




            idxs = [i for i in range(image.shape[0])]
            pool.map(slicepred, idxs)



            # for ii in range(image.shape[0]):
            #     imageslice = image[ii]
            #     shape_ = imageslice.shape
            #     # print(shape_)
                
            #     inputslice = normalize(imageslice)
            #     inputslice = cv2.resize(inputslice, image_size, cv2.INTER_NEAREST)
                
            #     transform = Compose([tf.ToTensor(),
            #             tf.Normalize(mean=[0.456], std=[0.224]),
            #             ])
            #     sample = {'image': inputslice}
            #     imagedata =transform(sample)
            #     imageslice = imagedata['image']
            #     predslice = self.step(imageslice, dataid)
            #     predslice = predslice.cpu().numpy()
            #     # plt.imshow(predslice)
            #     # plt.show()
            #     # plt.close()
            #     predslice = cv2.resize(predslice.astype(np.uint8), (shape_[1], shape_[0]), cv2.INTER_NEAREST)
                
            #     pred_mask[ii] = predslice

            #     if ii %10 == 0:
            #         print(f'curii {ii}')
            #         plt.imshow(image[ii])
            #         plt.title(f'{ii}')
            #         plt.show()
            #         plt.close()
            #         plt.imshow(pred_mask[ii])
            #         plt.show()
            #         plt.close()


            # pred_mask = tifffile.imread(f'F:\\salilab\\salilab\\projects\\auto-segmentation_Aneesh\\autoseg\\semantic\\results\\masks\\imageall_x_mito_predtiff\\GIP_30_1535.tiff')
            self.save_tiff(pred_mask, dataid)

            imagefile = dataid
            pred_list = pred_mask
            if self.type == 'imageall_x':
                #### We only have manual label in x view, so we evaluate here. Otherwise we save the tiff
                #### from three views to do postprocessing.
                self.gt = self.turnGtBigLabel(imagefile)
                dice, iou, acc = self.evaluate(pred_list, self.gt)
                # pred_list = pred_list.cpu().numpy()
                if args.mode == 'mito':
                    print("{}, mode: {}\n"
                        "mito: DICE: {:.2%} | IOU: {:.2%}".format(imagefile, self.type, dice[0], iou[0]))
                else:
                    print("{}, mode: {}\n"
                        "mem: DICE: {:.2%} | IOU: {:.2%}\n"
                        "nu: DICE: {:.2%} | IOU: {:.2%}".format(imagefile, self.type, dice[0], iou[0],dice[1], iou[1]))
                print('-------------------------------')



    def step(self, image, imagefile=None):
        self.imagefile= imagefile
        # print('line150', self.imagefile)
        image = torch.unsqueeze(image, 0)
        images = image
        # images = images.to(self.device)
        images = images.to(torch.device('cpu')) 
        # print('line167', images.shape)
        outputs = self.model(images)
        pred = torch.argmax(torch.softmax(outputs, dim=1), dim=1)
        if args.mode == 'mito':
            pred = self.turnPredMitoLabel(pred) #turn pred mito 0,1,2 to 0,1
        pred = torch.squeeze(pred, 0)
        # print('line161', images.shape)
        # temp_pred = pred.cpu().numpy()
        # plt.imshow(np.array(temp_img), cmap='gray')
        # plt.show()
        # plt.close()
        # plt.imshow(np.array(temp_pred))
        # plt.title(f'{idx}')
        # plt.show()
        # plt.close()



        # print('line163', pred_list)
        # print('line172', pred_list.shape)


        return pred


class Post_process(cellmapping_test):
    def __init__(self):
        super(Post_process, self).__init__(type='imageall_x')
        self.path = f'data/mask_xyz/'

    def getIsolateLabel(self, x, label):
        x_copy = x.copy()
        x_copy[x_copy != label] = 255
        x_copy[x_copy == label] = 1
        x_copy[x_copy == 255] = 0
        return x_copy

    def fix012(self, xyzLabel, xyz, x, y):
        if len(x[xyz==6] == 2) !=  len(x[xyz==6]) : #if 012 appear
            print('012 appear !')

    def save_tiff_coarse(self, pred, imagefile):
        # savepath = f'/group/xiangyi/iHuman-SIST/imageall_xyz_mask/fuse3D_mito_coarse_predtiff_0717fixflip'
        savepath = osp.join(f'{args.root_dir}//results//masks')

        make_dir(savepath)
        tifffile.imsave(f'{savepath}/{imagefile}_pred_mito.tiff', np.array(pred))
        print(f"save {savepath}/{imagefile}.tiff'")

    def cropCoarseLabel(self, label, imagefile):
        slice = label.shape[0]
        if glob.glob(osp.join(self.coarse_label_dir, imagefile, '*.png')) == []:
            print('did not have coarse label')
            return label
        mask = plt.imread(glob.glob(osp.join(self.coarse_label_dir, imagefile, '*.png'))[0], 0)
        maskrepeat = np.repeat(mask[None, ...], slice, axis=0)
        label[maskrepeat == 0] = 0
        return label

    def post_fusexyz(self):
        for imagefile in self.data_idxs:
            self.imagefile = imagefile
            print(self.imagefile)
            f'{args.root_dir}//results//masks//{self.type}_mito_predtiff'
            x = tifffile.imread(osp.join(f'{args.root_dir}//results//masks//imageall_x_{args.mode}_predtiff',   f'{self.imagefile}.tiff'))
            y = tifffile.imread(osp.join(f'{args.root_dir}//results//masks//imageall_y_{args.mode}_predtiff',   f'{self.imagefile}.tiff'))
            z = tifffile.imread(osp.join(f'{args.root_dir}//results//masks//imageall_z_{args.mode}_predtiff',   f'{self.imagefile}.tiff'))
            turned_y = np.rot90(y, k=1, axes=(2, 0))  ## y
            turned_z = np.rot90(z, k=3, axes=(1, 0))  ## z
            print('x y z turnedy turnedz', x.shape, y.shape,z.shape, turned_y.shape, turned_z.shape)
            assert x.shape == turned_y.shape and x.shape == turned_z.shape

            xyzLabel = np.zeros_like(x)
            for p in [0, 1, 2]:
                xyz = x + turned_y + turned_z
                xi = self.getIsolateLabel(x, p)  # x[x==p] = 1 x[x!=p] = 0
                yi = self.getIsolateLabel(turned_y, p)
                zi = self.getIsolateLabel(turned_z, p)
                xyzi = xi + yi + zi
                # only more than two location has this prediction,
                xyzLabel[xyzi >= 2] = p

            # If 012 case appears
            self.fix012(xyzLabel, xyz, x, turned_y)

            ############# Before Coarse Label #############
            if self.imagefile in self.test_idx:
                self.evaluate_again(xyzLabel, 'before')

            ## Optional: using the coarse label to refine the results. You can ignore this code if you don't want to use the extra labels.
            ## If you want to refine, the labels are in https://drive.google.com/drive/folders/12msPtKcQ7IPG7HD9OJZf5c5ZBntEiN1q?usp=sharing

            # self.cropCoarseLabel(xyzLabel, self.imagefile)
            # ############# After #############
            # if self.imagefile in self.test_idx:
            #     self.evaluate_again(xyzLabel, 'after')
            self.save_tiff_coarse(xyzLabel, self.imagefile)

    def evaluate_again(self, xyzLabel, flag='before'):
        gt = self.turnGtBigLabel(self.imagefile)
        dice, iou, _ = self.evaluate(xyzLabel, gt)
        if args.mode == 'mito':
            print("# {} coarse label refine \n"
                  "mito: DICE: {:.2%} | IOU: {:.2%}".format(flag, dice[0], iou[0]))
        else:
            print("# {} coarse label refine \n"
                  "mem: DICE: {:.2%} | IOU: {:.2%}\n"
                  "nu:  DICE: {:.2%} | IOU: {:.2%}".format(flag, dice[0], iou[0], dice[1], iou[1]))
        print('------------------------------- \n')


if __name__ == "__main__":
    pool = Pool(8)
    logging.basicConfig(level=logging.INFO)
    # for type in ['imageall_y', 'imageall_z']:
    for type in ['imageall_x', 'imageall_y', 'imageall_z']:
        print('#############', type, '###############')
        Test = cellmapping_test(type=type)
        Test.test()
    post_process = Post_process()
    post_process.post_fusexyz()


# %%

#%% test 



# %%

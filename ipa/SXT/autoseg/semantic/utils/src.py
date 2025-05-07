import os
import pickle
import numpy as np
from collections import OrderedDict
# from dice_loss import dice_coeff
import torch
import torch.nn.functional as F
from torch.autograd import Variable
import json
import ipdb
import matplotlib.pyplot as plt


def make_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def check_parallel(encoder_dict, decoder_dict):
    # check if the model was trained using multiple gpus
    trained_parallel = False
    for k, v in encoder_dict.items():
        if k[:7] == "module.":
            trained_parallel = True
        break
    if trained_parallel:
        # create new OrderedDict that does not contain "module."
        new_encoder_state_dict = OrderedDict()
        new_decoder_state_dict = OrderedDict()
        for k, v in encoder_dict.items():
            name = k[7:]  # remove "module."
            new_encoder_state_dict[name] = v
        for k, v in decoder_dict.items():
            name = k[7:]  # remove "module."
            new_decoder_state_dict[name] = v
        encoder_dict = new_encoder_state_dict
        decoder_dict = new_decoder_state_dict

    return encoder_dict, decoder_dict


def get_base_params(base_model, model):
    b = []
    if 'vgg' in base_model:
        b.append(model.base.features)
    else:
        b.append(model.base.conv1)
        b.append(model.base.bn1)
        b.append(model.base.layer1)
        b.append(model.base.layer2)
        b.append(model.base.layer3)
        b.append(model.base.layer4)

    for i in range(len(b)):
        for j in b[i].modules():
            jj = 0
            for k in j.parameters():
                jj += 1
                if k.requires_grad:
                    yield k


def get_skip_params(model):
    b = [model.sk1.parameters(), model.sk2.parameters(), model.sk3.parameters(), model.sk4.parameters(),
         model.sk5.parameters(), model.bn1.parameters(), model.bn2.parameters(), model.bn3.parameters(),
         model.bn4.parameters(), model.bn5.parameters()]

    for j in range(len(b)):
        for i in b[j]:
            yield i


def get_skip_dims(model_name):
    if model_name == 'resnet50' or model_name == 'resnet101':
        skip_dims_in = [2048, 1024, 512, 256, 64]
    elif model_name == 'resnet34':
        skip_dims_in = [512, 256, 128, 64, 64]
    elif model_name == 'vgg16':
        skip_dims_in = [512, 512, 256, 128, 64]

    return skip_dims_in


def center_crop(x, height=480, width=854):
    crop_h = torch.FloatTensor([x.size()[2]]).sub(height).div(-2)
    crop_w = torch.FloatTensor([x.size()[3]]).sub(width).div(-2)

    return F.pad(x, [
        crop_w.ceil().int()[0], crop_w.floor().int()[0],
        crop_h.ceil().int()[0], crop_h.floor().int()[0],
    ])


def get_optimizer(optim_name, lr, parameters, weight_decay=0, momentum=0.9):
    if optim_name == 'sgd':
        opt = torch.optim.SGD(filter(lambda p: p.requires_grad, parameters),
                              lr=lr, weight_decay=weight_decay, momentum=momentum)
    elif optim_name == 'adam':
        opt = torch.optim.Adam(filter(lambda p: p.requires_grad, parameters), lr=lr, weight_decay=weight_decay)
    elif optim_name == 'rmsprop':
        opt = torch.optim.RMSprop(filter(lambda p: p.requires_grad, parameters), lr=lr, weight_decay=weight_decay)
    return opt


def save_checkpoint(checkpoint_dir, encoder, decoder, enc_opt, dec_opt, paras):
    torch.save(encoder.state_dict(), os.path.join(checkpoint_dir, 'encoder.pt'))
    torch.save(decoder.state_dict(), os.path.join(checkpoint_dir, 'decoder.pt'))
    torch.save(enc_opt.state_dict(), os.path.join(checkpoint_dir, 'enc_opt.pt'))
    torch.save(dec_opt.state_dict(), os.path.join(checkpoint_dir, 'dec_opt.pt'))

    # save parameters for future use
    pickle.dump(paras, open(os.path.join(checkpoint_dir, 'paras.pkl'), 'wb'))


def save_checkpoint_lite(checkpoint_dir, model, optim, paras):
    torch.save(model.state_dict(), os.path.join(checkpoint_dir, 'model.pt'))
    torch.save(optim.state_dict(), os.path.join(checkpoint_dir, 'optim.pt'))
    # save parameters for future use
    pickle.dump(paras, open(os.path.join(checkpoint_dir, 'paras.pkl'), 'wb'))


def load_checkpoint(checkpoint_dir):
    encoder_dict = torch.load(os.path.join(checkpoint_dir, 'encoder.pt'))
    decoder_dict = torch.load(os.path.join(checkpoint_dir, 'decoder.pt'))
    enc_opt_dict = torch.load(os.path.join(checkpoint_dir, 'enc_opt.pt'))
    dec_opt_dict = torch.load(os.path.join(checkpoint_dir, 'dec_opt.pt'))

    # load parameters
    paras = pickle.load(open(os.path.join(checkpoint_dir, 'paras.pkl'), 'rb'))

    return encoder_dict, decoder_dict, enc_opt_dict, dec_opt_dict, paras


def load_checkpoint_lite(checkpoint_dir):
    state_dict = torch.load(os.path.join(checkpoint_dir, 'model.pt'))
    optim_dict = torch.load(os.path.join(checkpoint_dir, 'optim.pt'))

    return state_dict, optim_dict


def init_visdom(viz):
    mviz_pred = viz.image(np.zeros((480, 854)), opts=dict(title='Pred mask'))
    mviz_true = viz.image(np.zeros((480, 854)), opts=dict(title='True mask'))

    lot = viz.line(
        X=torch.zeros((1,)).cpu(),
        Y=torch.zeros((1,)).cpu(),
        opts=dict(
            xlabel='Iteration',
            ylabel='Loss',
            title='Running Loss',
            legend=['loss']
        )
    )

    elot = {}
    # epoch iou
    elot['iou'] = viz.line(
        X=torch.zeros((1,)).cpu(),
        Y=torch.zeros((1, 2)).cpu(),
        opts=dict(
            xlabel='Epoch',
            ylabel='IoU',
            title='IoU',
            legend=['train', 'val']
        )
    )

    # epoch loss
    elot['loss'] = viz.line(
        X=torch.zeros((1,)).cpu(),
        Y=torch.zeros((1, 2)).cpu(),
        opts=dict(
            xlabel='Epoch',
            ylabel='Loss',
            title='Total Loss',
            legend=['train', 'val']
        )
    )

    # text
    text = viz.text(text='start visdom')

    return lot, elot, mviz_pred, mviz_true, text


def IOU_dice(imPred, imLab, numClass):
    imPred = np.asarray(imPred.cpu()).copy()
    imLab = np.asarray(imLab.cpu()).copy()
    eps = 1
    # Remove classes from unlabeled pixels in gt image.
    # We should not penalize detections in unlabeled portions of the image.
    # imPred = imPred * (imLab > 0)

    # Compute area intersection:
    intersection = imPred * (imPred == imLab)
    (area_intersection, _) = np.histogram(intersection, bins=numClass - 1, range=(1, numClass - 1))

    # Compute area union:
    # We don't need background
    (area_pred, _) = np.histogram(imPred, bins=numClass - 1, range=(1, numClass - 1))  # 1, 2, 3, 4
    (area_lab, _) = np.histogram(imLab, bins=numClass - 1, range=(1, numClass - 1))
    area_union = area_pred + area_lab - area_intersection

    IOU = (area_intersection + eps) / (area_union + eps)
    Dice = (2 * area_intersection + eps) / (area_pred + area_lab + eps)

    return IOU, Dice


def dice_coeff_f(input, target):
    inter = torch.dot(input.view(-1), target.view(-1)) + 0.0001
    union = torch.sum(input) + torch.sum(target) + 0.0001

    t = 2 * inter.float() / union.float()
    return t


def accuracy(preds, label):
    valid = (label > 0)
    acc_sum = (valid * (preds == label)).sum()
    valid_sum = valid.sum()
    acc = float(acc_sum) / (float(valid_sum) + 1e-10)
    return acc


def evaluate_iou(x, y, num_class):
    """
    :param x: tensor (B, H, W) {0, 5} float
    :param y: tensor (B, 1, H, W) {0, 5} float
    :return: IOU: float
    """
    batch_size, h, w = x.size()

    x = x.view(batch_size, -1).long()
    y = y.view(batch_size, -1).long()
    acc = accuracy(x, y)
    IOU, Dice = IOU_dice(x, y, num_class)

    return IOU, acc, Dice


def class_balanced_cross_entropy_loss(outputs, label, size_average=True, batch_average=True):
    """Define the class balanced cross entropy loss to train the network
    Args:
    output: Output of the network
    label: Ground truth label
    Returns:
    Tensor that evaluates the loss
    """

    labels = label.float()

    batch_size = label.shape[0]

    for bid in range(batch_size):
        label = labels[bid]
        output = outputs[bid]
        label_map_bak = (label == 0).float()
        label_map_one = (label == 1).float()
        label_map_two = (label == 2).float()
        label_map_thr = (label == 3).float()
        label_map_for = (label == 4).float()
        num_back = torch.sum(label_map_bak)
        num_one = torch.sum(label == 1)
        num_ = torch.sum(label == 2)
        num_one = torch.sum(label == 3)
        num_one = torch.sum(label == 4)

    num_labels_pos = torch.sum(labels)
    num_labels_neg = torch.sum(1.0 - labels)
    num_total = num_labels_pos + num_labels_neg

    output_gt_zero = torch.ge(output, 0).float()
    loss_val = torch.mul(output, (labels - output_gt_zero)) - torch.log(
        1 + torch.exp(output - 2 * torch.mul(output, output_gt_zero)))

    loss_pos = torch.sum(-torch.mul(labels, loss_val))
    loss_neg = torch.sum(-torch.mul(1.0 - labels, loss_val))

    final_loss = num_labels_neg / num_total * loss_pos + num_labels_pos / num_total * loss_neg

    if size_average:
        final_loss /= int(np.prod(label.size()))
    elif batch_average:
        final_loss /= int(label.size()[0])

    return final_loss


def adjust_learning_rate(optimizer, epoch, step):
    """Sets the learning rate to the initial LR decayed by 10 every 30 epochs"""
    if np.mod(epoch + 1, step) == 0:
        for param_group in optimizer.param_groups:
            param_group['lr'] /= 10


def get_cityscapes_labels():
    return np.array([
        # [  0,   0,   0],
        [128, 64, 128],
        [244, 35, 232],
        [70, 70, 70],
        [102, 102, 156],
        [190, 153, 153],
        [153, 153, 153],
        [250, 170, 30],
        [220, 220, 0],
        [107, 142, 35],
        [152, 251, 152],
        [0, 130, 180],
        [220, 20, 60],
        [255, 0, 0],
        [0, 0, 142],
        [0, 0, 70],
        [0, 60, 100],
        [0, 80, 100],
        [0, 0, 230],
        [119, 11, 32]])


def get_pascal_labels():
    """Load the mapping that associates pascal classes with label colors
    Returns:
        np.ndarray with dimensions (21, 3)
    """
    return np.asarray([[0, 0, 0], [128, 0, 0], [0, 128, 0], [128, 128, 0],
                       [0, 0, 128], [128, 0, 128], [0, 128, 128], [128, 128, 128],
                       [64, 0, 0], [192, 0, 0], [64, 128, 0], [192, 128, 0],
                       [64, 0, 128], [192, 0, 128], [64, 128, 128], [192, 128, 128],
                       [0, 64, 0], [128, 64, 0], [0, 192, 0], [128, 192, 0],
                       [0, 64, 128]])


def encode_segmap(mask):
    """Encode segmentation label images as pascal classes
    Args:
        mask (np.ndarray): raw segmentation label image of dimension
          (M, N, 3), in which the Pascal classes are encoded as colours.
    Returns:
        (np.ndarray): class map with dimensions (M,N), where the value at
        a given location is the integer denoting the class index.
    """
    mask = mask.astype(int)
    label_mask = np.zeros((mask.shape[0], mask.shape[1]), dtype=np.int16)
    for ii, label in enumerate(get_pascal_labels()):
        label_mask[np.where(np.all(mask == label, axis=-1))[:2]] = ii
    label_mask = label_mask.astype(int)
    return label_mask


def decode_seg_map_sequence(label_masks, dataset='pascal'):
    rgb_masks = []
    for label_mask in label_masks:
        # print(label_mask.shape)
        rgb_mask = decode_segmap(label_mask, dataset)
        rgb_masks.append(rgb_mask)
    rgb_masks = torch.from_numpy(np.array(rgb_masks).transpose([0, 3, 1, 2]))
    return rgb_masks


def decode_segmap(label_mask, dataset, plot=False):
    """Decode segmentation class labels into a color image
    Args:
        label_mask (np.ndarray): an (M,N) array of integer values denoting
          the class label at each spatial location.
        plot (bool, optional): whether to show the resulting color image
          in a figure.
    Returns:
        (np.ndarray, optional): the resulting decoded color image.
    """
    if dataset == 'pascal':
        n_classes = 21
        label_colours = get_pascal_labels()
    elif dataset == 'cityscapes':
        n_classes = 19
        label_colours = get_cityscapes_labels()
    else:
        raise NotImplementedError

    r = label_mask.copy()
    g = label_mask.copy()
    b = label_mask.copy()
    # n_classes = np.unique(label_mask).shape[0]
    for ll in range(0, n_classes):
        r[label_mask == ll] = label_colours[ll, 0]
        g[label_mask == ll] = label_colours[ll, 1]
        b[label_mask == ll] = label_colours[ll, 2]
    rgb = np.zeros((label_mask.shape[0], label_mask.shape[1], 3))
    # print(rgb.shape, r.shape)
    rgb[:, :, 0] = r / 255.0
    rgb[:, :, 1] = g / 255.0
    rgb[:, :, 2] = b / 255.0
    if plot:
        plt.imshow(rgb)
        plt.show()
    else:
        return rgb
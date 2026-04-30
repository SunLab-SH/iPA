"""
Train SXT ISG Mask R-CNN Model

This script trains a Mask R-CNN model on the prepared COCO dataset.
"""

import os
import sys
import torch
import torchvision
from torchvision.models.detection import maskrcnn_resnet50_fpn
from torchvision.datasets import CocoDetection
from torchvision.transforms import functional as F
from pycocotools.coco import COCO
from pycocotools.cocoeval import COCOeval
from tqdm import tqdm
import numpy as np
from torch.utils.data import DataLoader

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

DATA_ROOT = 'data/sxt/sxt_isg_training/isg_maskrcnn_coco'
MODEL_OUTPUT_DIR = os.path.join(DATA_ROOT, 'models')
os.makedirs(MODEL_OUTPUT_DIR, exist_ok=True)

# Standard COCO transforms for Mask R-CNN
class Compose:
    def __init__(self, transforms):
        self.transforms = transforms

    def __call__(self, image, target):
        for t in self.transforms:
            image, target = t(image, target)
        return image, target

class ToTensor:
    def __call__(self, image, target):
        image = F.to_tensor(image)
        return image, target

def get_transform():
    return Compose([ToTensor()])

def get_model(num_classes=2):
    """Load a pre-trained Mask R-CNN and modify for our task."""
    model = maskrcnn_resnet50_fpn(weights='DEFAULT')
    # Replace the box predictor and mask predictor
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = torchvision.models.detection.faster_rcnn.FastRCNNPredictor(in_features, num_classes)
    
    in_features_mask = model.roi_heads.mask_predictor.conv5_mask.in_channels
    hidden_layer = 256
    model.roi_heads.mask_predictor = torchvision.models.detection.mask_rcnn.MaskRCNNPredictor(in_features_mask, hidden_layer, num_classes)
    
    # Disable internal resizing to match our 256x256 patches
    model.transform.min_size = (256,)
    model.transform.max_size = 256
    
    return model

def evaluate(model, data_loader, device, coco_gt_file):
    """Evaluate model on validation set using COCO mAP metrics."""
    model.eval()
    
    # Load ground truth
    coco_gt = COCO(coco_gt_file)
    
    # Collect all predictions
    coco_results = []
    
    with torch.no_grad():
        for images, targets in tqdm(data_loader, desc='Evaluating', leave=False):
            images = list(image.to(device) for image in images)
            predictions = model(images)
            
            for pred, target in zip(predictions, targets):
                # Get image_id - CocoDetection returns list of annotations
                img_id = None
                if isinstance(target, dict) and 'image_id' in target:
                    img_id = int(target['image_id'].item()) if hasattr(target['image_id'], 'item') else int(target['image_id'])
                elif isinstance(target, (list, tuple)) and len(target) > 0:
                    # Try to get image_id from first annotation
                    if isinstance(target[0], dict) and 'image_id' in target[0]:
                        img_id = int(target[0]['image_id'])
                
                if img_id is None:
                    continue
                
                boxes = pred['boxes'].cpu().numpy()
                scores = pred['scores'].cpu().numpy()
                labels = pred['labels'].cpu().numpy()
                
                for i in range(len(boxes)):
                    if scores[i] < 0.001:  # Lower threshold to capture early training predictions
                        continue
                    x1, y1, x2, y2 = boxes[i]
                    coco_results.append({
                        'image_id': img_id,
                        'category_id': int(labels[i]),
                        'bbox': [float(x1), float(y1), float(x2-x1), float(y2-y1)],
                        'score': float(scores[i])
                    })
    
    if len(coco_results) == 0:
        print("No predictions made")
        return 0.0
    
    # Save predictions to temporary file
    import json
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(coco_results, f)
        pred_file = f.name
    
    # Evaluate using COCO API
    try:
        coco_dt = coco_gt.loadRes(pred_file)
        coco_eval = COCOeval(coco_gt, coco_dt, 'bbox')
        coco_eval.evaluate()
        coco_eval.accumulate()
        coco_eval.summarize()
        
        # Get mAP @ IoU=0.50:0.95
        map_score = coco_eval.stats[0]
        print(f"\nmAP @ [0.50:0.95]: {map_score:.4f}")
        print(f"mAP @ 0.50: {coco_eval.stats[1]:.4f}")
        print(f"mAP @ 0.75: {coco_eval.stats[2]:.4f}")
        
        # Clean up
        import os
        os.unlink(pred_file)
        
        return map_score
    except Exception as e:
        print(f"Evaluation error: {e}")
        import os
        if os.path.exists(pred_file):
            os.unlink(pred_file)
        return 0.0

def train_one_epoch(model, optimizer, data_loader, device, epoch, scaler=None, print_freq=10):
    model.train()
    
    # Add tqdm progress bar
    pbar = tqdm(enumerate(data_loader), total=len(data_loader), desc=f'Epoch {epoch+1}')
    
    for i, (images, targets) in pbar:
        images = list(image.to(device) for image in images)
        
        # Convert COCO targets to Mask R-CNN format
        new_targets = []
        for t in targets:
            boxes = []
            labels = []
            masks = []
            
            # Get image size to create masks
            img_h, img_w = 256, 256
            
            for ann in t:
                x, y, w, h = ann['bbox']
                boxes.append([x, y, x + w, y + h])
                labels.append(1)
                
                # Create a simple binary mask from bbox (Optimized)
                mask = np.zeros((img_h, img_w), dtype=np.uint8)
                mask[int(y):int(y+h), int(x):int(x+w)] = 1
                masks.append(mask)
            
            if not boxes:
                boxes = [[0, 0, 1, 1]]
                labels = [1]
                masks = [np.zeros((img_h, img_w), dtype=np.uint8)]
            
            # Stack masks efficiently to avoid the warning
            masks_tensor = torch.as_tensor(np.stack(masks), dtype=torch.uint8).to(device)
                
            new_targets.append({
                "boxes": torch.as_tensor(boxes, dtype=torch.float32).to(device),
                "labels": torch.as_tensor(labels, dtype=torch.int64).to(device),
                "masks": masks_tensor,
            })
        targets = new_targets
        
        # Temporarily disable AMP to debug coordinate issues
        loss_dict = model(images, targets)
        losses = sum(loss for loss in loss_dict.values())
        
        optimizer.zero_grad()
        losses.backward()
        optimizer.step()
        
        # Update progress bar with current loss
        pbar.set_postfix({'loss': f'{losses.item():.4f}'})

def main(train_ratio=0.05, val_ratio=0.05):
    """
    Train Mask R-CNN on SXT ISG dataset.
    
    Args:
        train_ratio: Fraction of training data to use (default: 0.05 for quick debugging)
        val_ratio: Fraction of validation data to use (default: 0.05)
    """
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    print(f"Using device: {device}")
    print(f"Train ratio: {train_ratio*100:.0f}%, Val ratio: {val_ratio*100:.0f}%")
    
    # Use standard CocoDetection for training
    train_dataset_full = CocoDetection(
        root=os.path.join(DATA_ROOT, 'images'),
        annFile=os.path.join(DATA_ROOT, 'train.json'),
        transforms=get_transform()
    )
    
    # Subsample training data for debugging
    if train_ratio < 1.0:
        num_samples = int(len(train_dataset_full) * train_ratio)
        indices = list(range(num_samples))
        train_dataset = torch.utils.data.Subset(train_dataset_full, indices)
        print(f"Using {num_samples}/{len(train_dataset_full)} training samples ({train_ratio*100:.0f}%)")
    else:
        train_dataset = train_dataset_full
        print(f"Using all {len(train_dataset)} training samples")
    
    # Validation dataset with subsampling
    val_dataset_full = CocoDetection(
        root=os.path.join(DATA_ROOT, 'images'),
        annFile=os.path.join(DATA_ROOT, 'val.json'),
        transforms=get_transform()
    )
    
    if val_ratio < 1.0:
        num_val_samples = int(len(val_dataset_full) * val_ratio)
        val_indices = list(range(num_val_samples))
        val_dataset = torch.utils.data.Subset(val_dataset_full, val_indices)
        print(f"Using {num_val_samples}/{len(val_dataset_full)} validation samples ({val_ratio*100:.0f}%)")
    else:
        val_dataset = val_dataset_full
        print(f"Using all {len(val_dataset)} validation samples")
    
    train_loader = DataLoader(
        train_dataset, 
        batch_size=4,
        shuffle=True, 
        num_workers=4,
        collate_fn=lambda x: tuple(zip(*x))
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=4,
        shuffle=False,
        num_workers=4,
        collate_fn=lambda x: tuple(zip(*x))
    )
    
    model = get_model(num_classes=2)
    model.to(device)
    
    params = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.SGD(params, lr=0.005, momentum=0.9, weight_decay=0.0005)
    
    # Mixed Precision Training for Speed
    try:
        scaler = torch.cuda.amp.GradScaler()
    except:
        scaler = torch.amp.GradScaler('cuda')
    
    num_epochs = 20
    eval_interval = 2  # Evaluate every 2 epochs
    best_map = 0.0
    
    for epoch in range(num_epochs):
        print(f"\n=== Epoch {epoch+1}/{num_epochs} ===")
        train_one_epoch(model, optimizer, train_loader, device, epoch, scaler)
        
        # Evaluate on validation set every eval_interval epochs
        if (epoch + 1) % eval_interval == 0 or (epoch + 1) == num_epochs:
            val_map = evaluate(model, val_loader, device, os.path.join(DATA_ROOT, 'val.json'))
            print(f"Validation mAP: {val_map:.4f}")
            
            # Save best model
            if val_map > best_map:
                best_map = val_map
                torch.save(model.state_dict(), os.path.join(MODEL_OUTPUT_DIR, 'isg_mask_rcnn_best.pth'))
                print(f"Saved best model with mAP: {best_map:.4f}")
        else:
            print(f"Skipping evaluation (will evaluate at epoch {(epoch // eval_interval + 1) * eval_interval})")
        
        # Save checkpoint every 5 epochs
        if (epoch + 1) % 5 == 0:
            torch.save(model.state_dict(), os.path.join(MODEL_OUTPUT_DIR, f'isg_mask_rcnn_ep{epoch+1}.pth'))
        
    torch.save(model.state_dict(), os.path.join(MODEL_OUTPUT_DIR, 'isg_mask_rcnn_final.pth'))
    print(f"\nTraining complete! Best mAP: {best_map:.4f}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Train SXT ISG Mask R-CNN')
    parser.add_argument('--train-ratio', type=float, default=0.05,
                        help='Fraction of training data to use (default: 0.05)')
    parser.add_argument('--val-ratio', type=float, default=0.05,
                        help='Fraction of validation data to use (default: 0.05)')
    args = parser.parse_args()
    main(train_ratio=args.train_ratio, val_ratio=args.val_ratio)

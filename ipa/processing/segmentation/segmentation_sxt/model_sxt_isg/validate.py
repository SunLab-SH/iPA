import os
import sys
import time
import json
import warnings
from pathlib import Path
from datetime import datetime

import torch
import torch.nn.functional as F
import numpy as np
import tifffile
import mrcfile
from PIL import Image
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

from .unet import UNet
from .train import ISG3DDataset
from .utils.dice_score import dice_coeff, multiclass_dice_coeff

def calculate_metrics(y_true, y_pred, y_prob=None, class_names=None):
    """
    Calculate comprehensive evaluation metrics
    
    Args:
        y_true: Ground truth labels (numpy array)
        y_pred: Predicted labels (numpy array) 
        y_prob: Predicted probabilities (numpy array, optional)
        class_names: List of class names (optional)
    
    Returns:
        dict: Dictionary containing all metrics
    """
    metrics = {}
    
    # Basic classification metrics
    metrics['accuracy'] = accuracy_score(y_true, y_pred)
    metrics['precision'] = precision_score(y_true, y_pred, average='weighted', zero_division=0)
    metrics['recall'] = recall_score(y_true, y_pred, average='weighted', zero_division=0)
    metrics['f1_score'] = f1_score(y_true, y_pred, average='weighted', zero_division=0)
    
    # Per-class metrics
    if class_names is None:
        class_names = [f'Class_{i}' for i in range(len(np.unique(y_true)))]
    
    precision_per_class = precision_score(y_true, y_pred, average=None, zero_division=0)
    recall_per_class = recall_score(y_true, y_pred, average=None, zero_division=0)
    f1_per_class = f1_score(y_true, y_pred, average=None, zero_division=0)
    
    metrics['per_class_metrics'] = {}
    for i, class_name in enumerate(class_names):
        if i < len(precision_per_class):
            metrics['per_class_metrics'][class_name] = {
                'precision': float(precision_per_class[i]),
                'recall': float(recall_per_class[i]),
                'f1_score': float(f1_per_class[i])
            }
    
    # Confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    metrics['confusion_matrix'] = cm.tolist()
    
    # AUC metrics (if probabilities provided)
    if y_prob is not None:
        try:
            if len(np.unique(y_true)) == 2:  # Binary classification
                metrics['roc_auc'] = roc_auc_score(y_true, y_prob[:, 1] if y_prob.ndim > 1 else y_prob)
                
                # Precision-Recall AUC
                precision_curve, recall_curve, _ = precision_recall_curve(
                    y_true, y_prob[:, 1] if y_prob.ndim > 1 else y_prob
                )
                metrics['pr_auc'] = auc(recall_curve, precision_curve)
                
            else:  # Multi-class classification
                metrics['roc_auc'] = roc_auc_score(y_true, y_prob, multi_class='ovr', average='weighted')
                
        except Exception as e:
            print(f"Warning: Could not calculate AUC metrics: {e}")
            metrics['roc_auc'] = None
            metrics['pr_auc'] = None
    
    return metrics

def validate_model(model_path, 
                   test_img_dir, 
                   test_mask_dir, 
                   device='cuda',
                   scale_factor=0.5,
                   batch_size=8,
                   test_ids=None,
                   organelle='isg',
                   save_results=True,
                   results_dir=None):
    """
    Comprehensive model validation with detailed metrics
    
    Args:
        model_path: Path to trained model weights
        test_img_dir: Directory containing test images
        test_mask_dir: Directory containing test masks
        device: Device to run validation on
        scale_factor: Image scaling factor
        batch_size: Batch size for validation
        test_ids: List of test IDs to validate on
        organelle: Type of organelle being segmented
        save_results: Whether to save validation results
        results_dir: Directory to save results
    
    Returns:
        dict: Comprehensive validation results
    """
    
    print("=" * 60)
    print("STARTING MODEL VALIDATION")
    print("=" * 60)
    
    validation_start = time.time()
    
    # Setup results directory
    if results_dir is None:
        results_dir = Path(model_path).parent / 'validation_results'
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Load model
    print(f"Loading model from: {model_path}")
    
    # Determine number of classes from model file name or default
    if 'cell_nucleus' in str(model_path) or organelle == 'cell_nucleus':
        n_classes = 3
        class_names = ['background', 'cell', 'nucleus']
    else:
        n_classes = 2  # ISG, mito, etc.
        class_names = ['background', organelle]
    
    net = UNet(n_channels=3, n_classes=n_classes, bilinear=False)
    net.load_state_dict(torch.load(model_path, map_location=device))
    net.to(device)
    net.eval()
    
    print(f"Model loaded with {n_classes} classes: {class_names}")
    
    # Create test dataset
    print(f"Creating test dataset from:")
    print(f"  Images: {test_img_dir}")
    print(f"  Masks: {test_mask_dir}")
    
    test_dataset = ISG3DDataset(
        test_img_dir, 
        test_mask_dir, 
        scale=scale_factor,
        allowed_ids=test_ids
    )
    
    test_loader = DataLoader(
        test_dataset, 
        batch_size=batch_size, 
        shuffle=False, 
        num_workers=2, 
        pin_memory=True
    )
    
    print(f"Test dataset contains {len(test_dataset)} slices")
    
    # Validation loop
    all_predictions = []
    all_probabilities = []
    all_ground_truth = []
    all_dice_scores = []
    
    total_samples = 0
    
    print("Running validation...")
    
    with torch.no_grad():
        for batch_idx, batch in enumerate(test_loader):
            if (batch_idx + 1) % 10 == 0:
                print(f"Processing batch {batch_idx + 1}/{len(test_loader)}")
            
            images = batch['image'].to(device, dtype=torch.float32)
            true_masks = batch['mask'].to(device, dtype=torch.long)
            
            # Forward pass
            outputs = net(images)
            
            # Get probabilities and predictions
            if n_classes > 1:
                probs = F.softmax(outputs, dim=1)
                preds = torch.argmax(probs, dim=1)
            else:
                probs = torch.sigmoid(outputs)
                preds = (probs > 0.5).long()
            
            # Calculate Dice score for this batch
            if n_classes > 1:
                dice_score = multiclass_dice_coeff(probs, true_masks, reduce_batch_first=False)
            else:
                dice_score = dice_coeff(preds.float(), true_masks.float(), reduce_batch_first=False)
            
            all_dice_scores.append(dice_score.cpu().numpy())
            
            # Store predictions and ground truth
            batch_preds = preds.cpu().numpy().flatten()
            batch_true = true_masks.cpu().numpy().flatten()
            batch_probs = probs.cpu().numpy()
            
            all_predictions.extend(batch_preds)
            all_ground_truth.extend(batch_true)
            
            # Reshape probabilities for AUC calculation
            if n_classes > 1:
                # For multi-class, we need probabilities for each class
                batch_probs_reshaped = batch_probs.transpose(0, 2, 3, 1).reshape(-1, n_classes)
            else:
                # For binary, we just need the positive class probability
                batch_probs_reshaped = batch_probs.flatten()
            
            all_probabilities.extend(batch_probs_reshaped)
            
            total_samples += images.shape[0]
    
    # Convert to numpy arrays
    all_predictions = np.array(all_predictions)
    all_ground_truth = np.array(all_ground_truth)
    
    if n_classes > 1:
        all_probabilities = np.array(all_probabilities)
    else:
        all_probabilities = np.array(all_probabilities).reshape(-1, 1)
        # Add background probability for binary case
        all_probabilities = np.column_stack([1 - all_probabilities, all_probabilities])
    
    print(f"Validation completed on {total_samples} samples")
    
    # Calculate comprehensive metrics
    print("Calculating metrics...")
    metrics = calculate_metrics(
        all_ground_truth, 
        all_predictions, 
        all_probabilities,
        class_names
    )
    
    # Add Dice score metrics
    mean_dice = np.mean(all_dice_scores)
    std_dice = np.std(all_dice_scores)
    
    metrics['dice_score'] = {
        'mean': float(mean_dice),
        'std': float(std_dice),
        'min': float(np.min(all_dice_scores)),
        'max': float(np.max(all_dice_scores))
    }
    
    # Validation summary
    validation_time = time.time() - validation_start
    
    validation_results = {
        'model_path': str(model_path),
        'validation_time': validation_time,
        'total_samples': total_samples,
        'total_slices': len(test_dataset),
        'organelle': organelle,
        'n_classes': n_classes,
        'class_names': class_names,
        'test_ids': test_ids,
        'metrics': metrics,
        'validation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    # Print results summary
    print("\n" + "=" * 60)
    print("VALIDATION RESULTS SUMMARY")
    print("=" * 60)
    print(f"Model: {Path(model_path).name}")
    print(f"Organelle: {organelle}")
    print(f"Classes: {n_classes} ({', '.join(class_names)})")
    print(f"Validation time: {validation_time:.2f}s")
    print(f"Total samples: {total_samples}")
    print(f"Total slices: {len(test_dataset)}")
    print()
    print("OVERALL METRICS:")
    print(f"  Accuracy: {metrics['accuracy']:.4f}")
    print(f"  Precision (weighted): {metrics['precision']:.4f}")
    print(f"  Recall (weighted): {metrics['recall']:.4f}")
    print(f"  F1-Score (weighted): {metrics['f1_score']:.4f}")
    print(f"  Dice Score: {mean_dice:.4f} ± {std_dice:.4f}")
    
    if metrics.get('roc_auc') is not None:
        print(f"  ROC AUC: {metrics['roc_auc']:.4f}")
    if metrics.get('pr_auc') is not None:
        print(f"  PR AUC: {metrics['pr_auc']:.4f}")
    
    print("\nPER-CLASS METRICS:")
    for class_name, class_metrics in metrics['per_class_metrics'].items():
        print(f"  {class_name}:")
        print(f"    Precision: {class_metrics['precision']:.4f}")
        print(f"    Recall: {class_metrics['recall']:.4f}")
        print(f"    F1-Score: {class_metrics['f1_score']:.4f}")
    
    # Save results
    if save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_file = results_dir / f'validation_results_{organelle}_{timestamp}.json'
        
        with open(results_file, 'w') as f:
            json.dump(validation_results, f, indent=2)
        
        print(f"\nValidation results saved to: {results_file}")
        
        # Save confusion matrix plot
        plt.figure(figsize=(8, 6))
        cm = np.array(metrics['confusion_matrix'])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                   xticklabels=class_names, yticklabels=class_names)
        plt.title(f'Confusion Matrix - {organelle.upper()}')
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
        
        cm_file = results_dir / f'confusion_matrix_{organelle}_{timestamp}.png'
        plt.savefig(cm_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Confusion matrix saved to: {cm_file}")
    
    print("=" * 60)
    
    return validation_results

def validate_all_models(models_dir, 
                       test_img_dir, 
                       test_mask_dir,
                       test_ids=None,
                       results_dir=None):
    """
    Validate all models in a directory
    
    Args:
        models_dir: Directory containing model files
        test_img_dir: Directory containing test images
        test_mask_dir: Directory containing test masks
        test_ids: List of test IDs
        results_dir: Directory to save results
    """
    
    models_dir = Path(models_dir)
    model_files = list(models_dir.glob('*.pth'))
    
    print(f"Found {len(model_files)} model files in {models_dir}")
    
    all_results = {}
    
    for model_file in model_files:
        print(f"\n{'='*60}")
        print(f"Validating model: {model_file.name}")
        print(f"{'='*60}")
        
        # Determine organelle type from filename
        organelle = 'isg'  # default
        if 'mito' in model_file.name.lower():
            organelle = 'mito'
        elif 'cell' in model_file.name.lower() or 'nucleus' in model_file.name.lower():
            organelle = 'cell_nucleus'
        
        try:
            results = validate_model(
                model_path=model_file,
                test_img_dir=test_img_dir,
                test_mask_dir=test_mask_dir,
                test_ids=test_ids,
                organelle=organelle,
                results_dir=results_dir
            )
            
            all_results[model_file.name] = results
            
        except Exception as e:
            print(f"Error validating {model_file.name}: {e}")
            continue
    
    # Save combined results
    if results_dir and all_results:
        combined_file = Path(results_dir) / f'all_validation_results_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
        with open(combined_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        print(f"\nCombined results saved to: {combined_file}")
    
    return all_results

# if __name__ == "__main__":
#     # Example usage
#     import argparse
    
#     parser = argparse.ArgumentParser(description='Validate trained models')
#     parser.add_argument('--model', type=str, help='Path to model file')
#     parser.add_argument('--models_dir', type=str, help='Directory containing multiple models')
#     parser.add_argument('--test_img_dir', type=str, required=True, help='Test images directory')
#     parser.add_argument('--test_mask_dir', type=str, required=True, help='Test masks directory')
#     parser.add_argument('--results_dir', type=str, help='Results output directory')
#     parser.add_argument('--device', type=str, default='cuda', help='Device to use')
#     parser.add_argument('--batch_size', type=int, default=8, help='Batch size')
#     parser.add_argument('--organelle', type=str, default='isg', help='Organelle type')
    
#     args = parser.parse_args()
    
#     device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    
#     if args.model:
#         # Validate single model
#         validate_model(
#             model_path=args.model,
#             test_img_dir=args.test_img_dir,
#             test_mask_dir=args.test_mask_dir,
#             device=device,
#             batch_size=args.batch_size,
#             organelle=args.organelle,
#             results_dir=args.results_dir
#         )
#     elif args.models_dir:
#         # Validate all models in directory
#         validate_all_models(
#             models_dir=args.models_dir,
#             test_img_dir=args.test_img_dir,
#             test_mask_dir=args.test_mask_dir,
#             results_dir=args.results_dir
#         )
#     else:
#         print("Please specify either --model or --models_dir")

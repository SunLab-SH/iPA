import os
import numpy as np
import torch
import torch.nn.functional as F
from PIL import Image
from torchvision import transforms
import pandas as pd
import tifffile
import mrcfile
import warnings
from .unet import UNet
from .utils.data_loading import BasicDataset

# Filter MRC file warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="mrcfile")
warnings.filterwarnings("ignore", message=".*Unrecognised machine stamp.*")
warnings.filterwarnings("ignore", message=".*Map ID string not found.*")
warnings.filterwarnings("ignore", message=".*data block cannot be read.*")

filepath = os.path.dirname(os.path.abspath(__file__))
# print('filepath:', filepath)

def isg_mask_pred(net, img3d, device, scale_factor=0.5, out_threshold=0.5):
    net.eval()
    outputmask = np.zeros_like(img3d)
    
    total_frames = img3d.shape[0]
    for i in range(total_frames):
        # Print progress for every frame
        if (i + 1) % 50 == 0 or i == total_frames - 1:
            print(f"Processing frame {i+1}/{total_frames}")
        
        img_slice = img3d[i]
        img_rgb = np.stack([img_slice]*3, axis=-1).astype(np.uint8)
        pil_img = Image.fromarray(img_rgb)

        img_tensor = torch.from_numpy(BasicDataset.preprocess(pil_img, scale_factor, is_mask=False)).unsqueeze(0)
        img_tensor = img_tensor.to(device=device, dtype=torch.float32)

        with torch.no_grad():
            output = net(img_tensor)
            if net.n_classes > 1:
                probs = F.softmax(output, dim=1)[0]
            else:
                probs = torch.sigmoid(output)[0]

            tf = transforms.Compose([
                transforms.ToPILImage(),
                transforms.Resize((pil_img.size[1], pil_img.size[0])),
                transforms.ToTensor()
            ])
            full_mask = tf(probs.cpu()).squeeze()

            if net.n_classes == 1:
                mask_np = (full_mask > out_threshold).numpy()
            else:
                mask_np = F.one_hot(full_mask.argmax(dim=0), net.n_classes).permute(2, 0, 1).numpy()
            
            if net.n_classes == 1:
                outputmask[i] = mask_np
            else:
                outputmask[i] = mask_np[1]

    return outputmask


def load_lac_factor():
    # lac_path = r'./Autosegmentation LIst.xlsx'
    # df = pd.read_excel(lac_path)
    # for _, row in df.iterrows():
    #     if idx in str(row['Cell ID']):
    #         return row['LAC Normalizing Factor']
    return 33.33


def convert_img(img_3d):
    min_val = np.min(img_3d)
    max_val = 1  # 固定值
    img_3d_new = (img_3d - min_val) / max_val * 255
    return img_3d_new.astype(np.uint8)


def load_img(path):

    img_path = path
    if img_path.endswith('.tif'):
        img = tifffile.imread(img_path)
    else:
        img = mrcfile.open(img_path, permissive=True).data

    lac_factor = load_lac_factor()
    img = img * lac_factor
    img[img > 1] = 1

    img = convert_img(img)
    return img


def run_sphere_like_organelle_segmentation(save_dir=None, 
                                         mode='isg',
                                         pool_processes=6,
                                         dataid=None,
                                         image_data=None,
                                         scale_factor=0.5,
                                         out_threshold=0.5,
                                         lac_factor=33.33,
                                         model_path=None,
                                         ground_truth_data=None,  # New parameter for validation
                                         calculate_metrics=False):  # New parameter
    """
    Run sphere-like organelle segmentation with batch processing support
    
    Args:
        save_dir: str, directory to save results
        mode: str, segmentation mode (default: 'isg')
        pool_processes: int, number of processes for parallel processing
        dataid: list, list of data IDs to process
        image_data: dict, dictionary mapping dataid to image data or file paths
        scale_factor: float, scale factor for model input (default: 0.5)
        out_threshold: float, threshold for binary classification (default: 0.5)
        lac_factor: float, LAC normalizing factor (default: 33.33)
        model_path: str, path to model weights (optional, uses default if None)
        ground_truth_data: dict, dictionary mapping dataid to ground truth masks (optional)
        calculate_metrics: bool, whether to calculate validation metrics (default: False)
    
    Returns:
        dict: Results including predictions and metrics (if calculated)
    """
    
    # Set default model path if not provided
    if model_path is None:
        model_path = f'{filepath}/vesicle_mask0.63_processed.pth'
    
    # Set default save directory if not provided
    if save_dir is None:
        save_dir = f'{filepath}/results/isg_semantic_mask'
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load model
    net = UNet(n_channels=3, n_classes=2, bilinear=False)
    net.load_state_dict(torch.load(model_path, map_location=device))
    net.to(device)
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Handle single data processing (backward compatibility)
    if dataid is None or image_data is None:
        print("Warning: Using legacy single data processing mode")
        return {}
    
    results = {
        'predictions': {},
        'metrics': {} if calculate_metrics else None
    }
    
    # Batch processing
    if isinstance(dataid, list):
        from multiprocessing import Pool
        
        def process_single_data(data_id):
            try:
                # Get image data
                if isinstance(image_data[data_id], str):
                    # If path string, load image
                    img3d = load_img_with_params(image_data[data_id], lac_factor)
                else:
                    # If already loaded image data
                    img3d = convert_img_with_params(image_data[data_id], lac_factor)
                
                # Run prediction
                mask = isg_mask_pred(net, img3d, device, scale_factor, out_threshold)
                
                # Save result
                output_file = os.path.join(save_dir, f'{data_id}_{mode}_mask.tif')
                tifffile.imsave(output_file, mask)
                print(f"Saved prediction mask for {data_id} to {output_file}")
                
                # Calculate metrics if ground truth provided
                metrics = None
                if calculate_metrics and ground_truth_data and data_id in ground_truth_data:
                    metrics = calculate_prediction_metrics(mask, ground_truth_data[data_id])
                    print(f"Metrics for {data_id}: Dice={metrics['dice_score']:.4f}, IoU={metrics['iou']:.4f}")
                
                return data_id, True, mask, metrics
                
            except Exception as e:
                print(f"Error processing {data_id}: {str(e)}")
                return data_id, False, None, None
        
        # Use multiprocessing for batch processing
        if pool_processes > 1:
            with Pool(processes=pool_processes) as pool:
                process_results = pool.map(process_single_data, dataid)
        else:
            process_results = [process_single_data(data_id) for data_id in dataid]
        
        # Collect results
        successful = 0
        for data_id, success, mask, metrics in process_results:
            if success:
                successful += 1
                results['predictions'][data_id] = mask
                if metrics and calculate_metrics:
                    results['metrics'][data_id] = metrics
        
        print(f"Successfully processed {successful}/{len(dataid)} datasets")
        
        # Calculate aggregate metrics
        if calculate_metrics and results['metrics']:
            aggregate_metrics = calculate_aggregate_metrics(results['metrics'])
            results['aggregate_metrics'] = aggregate_metrics
            print(f"Aggregate metrics: Dice={aggregate_metrics['mean_dice']:.4f}±{aggregate_metrics['std_dice']:.4f}")
        
    else:
        # Single data processing
        data_id = dataid
        if isinstance(image_data[data_id], str):
            img3d = load_img_with_params(image_data[data_id], lac_factor)
        else:
            img3d = convert_img_with_params(image_data[data_id], lac_factor)
        
        mask = isg_mask_pred(net, img3d, device, scale_factor, out_threshold)
        
        output_file = os.path.join(save_dir, f'{data_id}_{mode}_mask.tif')
        tifffile.imsave(output_file, mask)
        print(f"Saved prediction mask for {data_id} to {output_file}")
        
        results['predictions'][data_id] = mask
        
        # Calculate metrics if ground truth provided
        if calculate_metrics and ground_truth_data and data_id in ground_truth_data:
            metrics = calculate_prediction_metrics(mask, ground_truth_data[data_id])
            results['metrics'] = {data_id: metrics}
            print(f"Metrics for {data_id}: Dice={metrics['dice_score']:.4f}, IoU={metrics['iou']:.4f}")
    
    return results


def calculate_prediction_metrics(pred_mask, gt_mask):
    """
    Calculate metrics between predicted and ground truth masks
    
    Args:
        pred_mask: numpy array, predicted mask
        gt_mask: numpy array, ground truth mask
    
    Returns:
        dict: Dictionary containing calculated metrics
    """
    import numpy as np
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
    
    # Ensure masks are binary
    pred_binary = (pred_mask > 0).astype(np.uint8)
    gt_binary = (gt_mask > 0).astype(np.uint8)
    
    # Flatten for sklearn metrics
    pred_flat = pred_binary.flatten()
    gt_flat = gt_binary.flatten()
    
    # Calculate metrics
    metrics = {}
    
    # Basic metrics
    metrics['accuracy'] = float(accuracy_score(gt_flat, pred_flat))
    metrics['precision'] = float(precision_score(gt_flat, pred_flat, zero_division=0))
    metrics['recall'] = float(recall_score(gt_flat, pred_flat, zero_division=0))
    metrics['f1_score'] = float(f1_score(gt_flat, pred_flat, zero_division=0))
    
    # Dice score
    intersection = np.sum(pred_binary * gt_binary)
    dice_score = (2.0 * intersection) / (np.sum(pred_binary) + np.sum(gt_binary) + 1e-8)
    metrics['dice_score'] = float(dice_score)
    
    # IoU (Jaccard index)
    union = np.sum((pred_binary + gt_binary) > 0)
    iou = intersection / (union + 1e-8)
    metrics['iou'] = float(iou)
    
    # Sensitivity (True Positive Rate)
    tp = np.sum(pred_binary * gt_binary)
    fn = np.sum(gt_binary * (1 - pred_binary))
    sensitivity = tp / (tp + fn + 1e-8)
    metrics['sensitivity'] = float(sensitivity)
    
    # Specificity (True Negative Rate)
    tn = np.sum((1 - pred_binary) * (1 - gt_binary))
    fp = np.sum(pred_binary * (1 - gt_binary))
    specificity = tn / (tn + fp + 1e-8)
    metrics['specificity'] = float(specificity)
    
    return metrics


def calculate_aggregate_metrics(metrics_dict):
    """
    Calculate aggregate metrics across multiple samples
    
    Args:
        metrics_dict: dict, dictionary of metrics for each sample
    
    Returns:
        dict: Aggregate statistics
    """
    import numpy as np
    
    if not metrics_dict:
        return {}
    
    # Collect all metric values
    metric_names = list(next(iter(metrics_dict.values())).keys())
    aggregated = {}
    
    for metric_name in metric_names:
        values = [metrics[metric_name] for metrics in metrics_dict.values()]
        values = np.array(values)
        
        aggregated[f'mean_{metric_name}'] = float(np.mean(values))
        aggregated[f'std_{metric_name}'] = float(np.std(values))
        aggregated[f'min_{metric_name}'] = float(np.min(values))
        aggregated[f'max_{metric_name}'] = float(np.max(values))
        aggregated[f'median_{metric_name}'] = float(np.median(values))
    
    return aggregated

# Keep original function for backward compatibility
def run_sphere_like_organelle_segmentation_legacy(datapath=None, cell_id=None, output_dir=None):
    """Legacy function for backward compatibility"""
    model_path = f'{filepath}/vesicle_mask0.63_processed.pth'
    if output_dir is None:
        output_dir = f'{filepath}/results/isg_semantic_mask'
    
    # Convert to new format
    run_sphere_like_organelle_segmentation(
        save_dir=output_dir,
        mode='isg',
        pool_processes=1,
        dataid=[cell_id],
        image_data={cell_id: datapath}
    )

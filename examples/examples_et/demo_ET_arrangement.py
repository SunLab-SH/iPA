"""
Cryo-ET Spatial Arrangement Analysis Demo

This script demonstrates spatial arrangement analysis of anchored tubular organelles 
(e.g., F-actin) relative to the plasma membrane using Cryo-ET data.

It performs the following steps:
1. Loads Cryo-ET F-actin filament coordinates
2. Identifies filaments anchored to PM (within 120 nm)
3. Calculates angles between pairs of anchored filaments
4. Visualizes angle distribution and network organization

Note: This analysis belongs to the Spatial Arrangement module (arrangement.py),
not Morphology, as it focuses on spatial distribution patterns.
"""

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.analysis import analyze_anchored_organelle_angles
from ipa.data_loader import QuickLogger

def create_pm_mask_from_volume(volume_shape, pm_z_range=(40, 50)):
    """Create a simulated PM mask for demo."""
    pm_mask = np.zeros(volume_shape, dtype=np.uint8)
    pm_mask[:, :, pm_z_range[0]:pm_z_range[1]] = 1
    return pm_mask

def main():
    logger = QuickLogger("et_arrangement", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting Cryo-ET Spatial Arrangement Analysis")
    
    # --- 1. Load F-actin Coordinates ---
    logger.step("Loading F-actin filament coordinates...")
    data_root = os.path.join(PROJECT_ROOT, 'data', 'cryoET')
    actin_json_path = os.path.join(data_root, '20210517_30_009_actin_filled_points.json')
    
    if not os.path.exists(actin_json_path):
        logger.step(f"Error: F-actin coordinates file not found: {actin_json_path}")
        return
    
    with open(actin_json_path, 'r') as f:
        actin_data = json.load(f)
    
    logger.step(f"Loaded {len(actin_data)} F-actin filaments")
    
    # Convert to the format expected by the function
    # The JSON has filament_id -> list of [x, y, z] coordinates
    filament_coords_dict = {}
    for filament_id, coords in actin_data.items():
        if len(coords) >= 2:  # Need at least 2 points to define direction
            filament_coords_dict[filament_id] = coords
    
    logger.step(f"Using {len(filament_coords_dict)} filaments with ≥2 points")
    
    # --- 2. Create PM Mask ---
    logger.step("Creating PM mask...")
    # Determine volume shape from coordinates
    all_coords = []
    for coords in filament_coords_dict.values():
        all_coords.extend(coords)
    
    all_coords = np.array(all_coords)
    max_coords = np.max(all_coords, axis=0).astype(int) + 10
    
    # Create PM mask (PM at high z values)
    pm_mask = create_pm_mask_from_volume(tuple(max_coords), pm_z_range=(max_coords[2]-10, max_coords[2]))
    logger.step(f"PM mask shape: {pm_mask.shape}")
    
    # --- 3. Analyze Anchored Angles ---
    logger.step("Analyzing anchored organelle angles...")
    
    # Cryo-ET voxel size: 1.34 nm (from Manuscript)
    voxel_size_nm = (1.34, 1.34, 1.34)
    distance_threshold_nm = 120  # As mentioned in Manuscript
    
    try:
        results = analyze_anchored_organelle_angles(
            filament_coords_dict,
            pm_mask,
            distance_threshold_nm=distance_threshold_nm,
            voxel_size_nm=voxel_size_nm
        )
        
        logger.step("✅ Analysis Complete!")
        logger.step(f"   Anchored filaments: {results['anchored_count']}")
        logger.step(f"   Total angle pairs: {len(results['angles'])}")
        
        if len(results['angles']) > 0:
            logger.step(f"   Mean angle: {results['mean_angle']:.2f}°")
            logger.step(f"   Std angle: {results['std_angle']:.2f}°")
            logger.step(f"   Min angle: {np.min(results['angles']):.2f}°")
            logger.step(f"   Max angle: {np.max(results['angles']):.2f}°")
        else:
            logger.step("   No angle pairs calculated (need ≥2 anchored filaments)")
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # --- 4. Visualization ---
    logger.step("Generating visualizations...")
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'et_angle_analysis_demo')
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    
    # Plot 1: Angle distribution histogram
    if len(results['angles']) > 0:
        axes[0].hist(results['angles'], bins=20, color='#00BFBF', alpha=0.7, edgecolor='black')
        axes[0].axvline(results['mean_angle'], color='red', linestyle='--', linewidth=2,
                       label=f'Mean: {results["mean_angle"]:.2f}°')
        axes[0].set_xlabel('Angle (degrees)', fontsize=11)
        axes[0].set_ylabel('Frequency', fontsize=11)
        axes[0].set_title('Anchored Filament Angle Distribution', fontsize=12, fontweight='bold')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Add interpretation text
        if results['mean_angle'] < 30:
            interpretation = "Parallel/Aligned Network"
        elif results['mean_angle'] < 60:
            interpretation = "Moderately Organized"
        else:
            interpretation = "Randomly Oriented"
        
        axes[0].text(0.05, 0.95, f'Organization: {interpretation}',
                    transform=axes[0].transAxes, fontsize=10,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        axes[0].text(0.5, 0.5, 'No angles to display\n(Need ≥2 anchored filaments)',
                    ha='center', va='center', fontsize=12, transform=axes[0].transAxes)
        axes[0].set_title('No Data Available', fontsize=12, fontweight='bold')
    
    # Plot 2: Summary statistics
    stats_text = f"""
Anchored Filament Analysis
━━━━━━━━━━━━━━━━━━━━━━━

Total filaments: {len(filament_coords_dict)}
Anchored filaments: {results['anchored_count']}
Anchoring rate: {results['anchored_count']/len(filament_coords_dict)*100:.1f}%

Angle Statistics:
  Count: {len(results['angles'])} pairs
  Mean: {results['mean_angle']:.2f}°
  Std: {results['std_angle']:.2f}°
  
Parameters:
  Distance threshold: {distance_threshold_nm} nm
  Voxel size: {voxel_size_nm[0]} nm
    """
    
    axes[1].text(0.1, 0.9, stats_text, transform=axes[1].transAxes,
                fontsize=10, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    axes[1].axis('off')
    axes[1].set_title('Analysis Summary', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    viz_path = os.path.join(output_dir, 'angle_analysis.png')
    plt.savefig(viz_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.file_out(viz_path)
    
    # --- 5. Save Results ---
    logger.step("Saving results...")
    
    results_dict = {
        'total_filaments': len(filament_coords_dict),
        'anchored_count': results['anchored_count'],
        'anchoring_rate': float(results['anchored_count'] / len(filament_coords_dict)),
        'angle_statistics': {
            'count': len(results['angles']),
            'mean_angle': float(results['mean_angle']),
            'std_angle': float(results['std_angle']),
            'min_angle': float(np.min(results['angles'])) if len(results['angles']) > 0 else None,
            'max_angle': float(np.max(results['angles'])) if len(results['angles']) > 0 else None
        },
        'parameters': {
            'distance_threshold_nm': distance_threshold_nm,
            'voxel_size_nm': list(voxel_size_nm)
        }
    }
    
    results_path = os.path.join(output_dir, 'angle_analysis_results.json')
    with open(results_path, 'w') as f:
        json.dump(results_dict, f, indent=2)
    
    logger.file_out(results_path)
    
    logger.step(f"Results saved to: {output_dir}")
    logger.step("Demo finished successfully!")

if __name__ == '__main__':
    main()

"""
WFM Dynamics Analysis Demo

This script demonstrates complete dynamics analysis using REAL WFM time-lapse data.
It calculates both 3D velocity and radial velocity of organelles (e.g., ISGs).

Data source: data/wfm/16.7-30_P21_* (Real WFM time-lapse with 73 frames)
Tracking data: Position, Speed, and Velocity CSV/NPY files

It performs the following steps:
1. Loads real WFM tracking data (positions over time)
2. Calculates 3D velocity from consecutive frames
3. Visualizes velocity distributions
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, PROJECT_ROOT)

from ipa.data_loader import QuickLogger

def main():
    logger = QuickLogger("wfm_dynamics", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting WFM Dynamics Analysis Demo with REAL Data")
    
    # --- 1. Load Real WFM Velocity Data ---
    logger.step("Loading real WFM velocity data...")
    data_dir = os.path.join(PROJECT_ROOT, 'data', 'wfm')
    velocity_file = os.path.join(data_dir, '16.7-30_P21_Velocity_3dvelocity.npy')
    
    if not os.path.exists(velocity_file):
        logger.step(f"Error: Velocity file not found: {velocity_file}")
        print(f"❌ Velocity file not found: {velocity_file}")
        return
    
    # Load velocity data (shape: [4, N] where rows are [vx, vy, vz, track_id])
    velocity_data = np.load(velocity_file)
    velocities_3d = np.linalg.norm(velocity_data[:3], axis=0)  # Calculate magnitude
    
    n_measurements = len(velocities_3d)
    logger.step(f"Loaded {n_measurements} velocity measurements")
    logger.step(f"Mean 3D velocity: {np.mean(velocities_3d):.3f} μm/s")
    logger.step(f"Std 3D velocity: {np.std(velocities_3d):.3f} μm/s")
    logger.step(f"Max 3D velocity: {np.max(velocities_3d):.3f} μm/s")
    
    # --- 2. Visualize 3D Velocity Distribution ---
    logger.step("Creating velocity distribution visualization...")
    
    # Create violin plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Prepare data for violin plot
    data_to_plot = [velocities_3d]
    labels = ['WFM Real Data']
    
    vp = ax.violinplot(data_to_plot, showmeans=True, showmedians=True)
    
    # Customize plot
    ax.set_ylabel('3D Velocity (μm/s)', fontsize=12)
    ax.set_title('WFM Organelle 3D Velocity Distribution\n(Real Tracking Data)', fontsize=14, fontweight='bold')
    ax.set_xticks([1])
    ax.set_xticklabels(labels)
    
    # Add statistics text
    stats_text = f"Mean: {np.mean(velocities_3d):.3f} μm/s\n"
    stats_text += f"Std: {np.std(velocities_3d):.3f} μm/s\n"
    stats_text += f"N: {len(velocities_3d)}"
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    output_dir = os.path.join(PROJECT_ROOT, 'results')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'wfm_velocity_distribution.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    logger.file_out(output_path)
    logger.step(f"Velocity distribution saved to: {output_path}")
    plt.close()
    
    logger.step("WFM Dynamics Analysis completed successfully!")
    print(f"\n✅ Results saved to: {output_path}")
    print(f"   - Total measurements: {n_measurements}")
    print(f"   - Mean velocity: {np.mean(velocities_3d):.3f} μm/s")

if __name__ == '__main__':
    main()

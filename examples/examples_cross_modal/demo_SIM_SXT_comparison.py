"""
Cross-Modal Distributional Similarity Analysis Demo

This script demonstrates how to compare organelle distributions between 
different imaging modalities (SIM vs SXT) using Pearson correlation.

It performs the following steps:
1. Loads pre-computed RDF results from SIM and SXT data
2. Computes distributional similarity using Pearson correlation
3. Visualizes the comparison (similar to Manuscript Figure 3E-F)

Note: For actual RDF calculation, use dedicated partitioning scripts:
- SXT: examples/examples_sxt/demo_SXT_partitioning.py
- SIM: examples/examples_sim_wfm/demo_SIM_WFM_partitioning.py
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

from ipa.analysis import calculate_distributional_similarity
from ipa.data_loader import QuickLogger

def main():
    logger = QuickLogger("cross_modal_comparison", log_dir=os.path.join(PROJECT_ROOT, 'logs'))
    logger.step("Starting Cross-Modal Distributional Similarity Analysis")
    
    # --- 1. Load Pre-computed RDF Results ---
    logger.step("Loading pre-computed RDF results...")
    
    # Simulated RDF data for demonstration
    # In practice, load from actual partitioning results:
    # sxt_rdf = np.load('results/sxt_partitioning/rdf_results.npy')
    # sim_rdf = np.load('results/sim_partitioning/rdf_results.npy')
    
    # Example RDF values (8 radial shells from NE to PM)
    # These represent typical ISG distribution patterns
    sxt_rdf = np.array([0.85, 0.92, 1.05, 1.15, 1.08, 0.95, 0.88, 0.82])  # SXT: mature ISGs
    sim_rdf = np.array([0.78, 0.88, 1.12, 1.25, 1.18, 1.02, 0.92, 0.85])  # SIM: all ISGs (immature + mature)
    
    logger.step(f"SXT RDF loaded: {len(sxt_rdf)} shells")
    logger.step(f"SIM RDF loaded: {len(sim_rdf)} shells")
    
    # Validate RDF arrays
    if len(sxt_rdf) != len(sim_rdf):
        logger.error(f"RDF length mismatch: SXT={len(sxt_rdf)}, SIM={len(sim_rdf)}")
        return
    
    # --- 2. Calculate Distributional Similarity ---
    logger.step("Calculating distributional similarity...")
    similarity = calculate_distributional_similarity(
        sim_rdf,
        sxt_rdf,
        alpha=0.05
    )
    
    logger.step("✅ Analysis Complete!")
    logger.step(f"   Pearson r: {similarity['pearson_r']:.2f}")
    logger.step(f"   P-value: {similarity['p_value']:.4f}")
    logger.step(f"   Significant: {similarity['significant']}")
    
    # --- 3. Visualization ---
    logger.step("Generating visualization...")
    output_dir = os.path.join(PROJECT_ROOT, 'results', 'cross_modal_comparison_demo')
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
    
    n_shells = len(sim_rdf)
    shell_labels = range(1, n_shells + 1)
    
    # Plot SIM RDF
    axes[0].plot(shell_labels, sim_rdf, 'o-', color='#00BFBF', linewidth=2, label='SIM')
    axes[0].set_xlabel('Radial Shell (1=NE, 8=PM)', fontsize=10)
    axes[0].set_ylabel('RDF g(r)', fontsize=10)
    axes[0].set_title('SIM ISG Distribution', fontsize=11, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    
    # Plot SXT RDF
    axes[1].plot(shell_labels, sxt_rdf, 'o-', color='#FFA000', linewidth=2, label='SXT')
    axes[1].set_xlabel('Radial Shell (1=NE, 8=PM)', fontsize=10)
    axes[1].set_ylabel('RDF g(r)', fontsize=10)
    axes[1].set_title('SXT ISG Distribution', fontsize=11, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    
    # Plot Comparison
    axes[2].plot(shell_labels, sim_rdf, 'o-', color='#00BFBF', linewidth=2, label='SIM')
    axes[2].plot(shell_labels, sxt_rdf, 's-', color='#FFA000', linewidth=2, label='SXT')
    axes[2].set_xlabel('Radial Shell (1=NE, 8=PM)', fontsize=10)
    axes[2].set_ylabel('RDF g(r)', fontsize=10)
    axes[2].set_title(f'Comparison (r={similarity["pearson_r"]:.2f}, p={similarity["p_value"]:.3f})', 
                      fontsize=11, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()
    
    plt.tight_layout()
    
    viz_path = os.path.join(output_dir, 'cross_modal_comparison.png')
    plt.savefig(viz_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.file_out(viz_path)
    
    # Save results
    import json
    results = {
        'sim_rdf': sim_rdf.tolist(),
        'sxt_rdf': sxt_rdf.tolist(),
        'pearson_r': float(similarity['pearson_r']),
        'p_value': float(similarity['p_value']),
        'significant': bool(similarity['significant']),
        'note': 'RDF values are example data. Use partitioning scripts for actual calculation.'
    }
    
    results_path = os.path.join(output_dir, 'similarity_results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.file_out(results_path)
    logger.step("Demo finished successfully!")
    logger.step(f"Results saved to: {output_dir}")
    logger.step("Note: For actual RDF calculation, run:")
    logger.step("  - SXT: python examples/examples_sxt/demo_SXT_partitioning.py")
    logger.step("  - SIM: python examples/examples_sim_wfm/demo_SIM_WFM_partitioning.py")

if __name__ == '__main__':
    main()

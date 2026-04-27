import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
import plotly



def visualize_partitions(partition_masks, slice_idx=None, save_path=None):
    """
    Visualize radial partition zones in 2D or 3D.
    
    Parameters
    ----------
    partition_masks : list of numpy.ndarray
        List of partition masks for each radial zone
    slice_idx : int, optional
        Z-slice index for 2D visualization (shows 3D if None)
    save_path : str, optional
        Path to save the visualization
    """
    if slice_idx is not None:
        # 2D visualization
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create composite image with different colors for each zone
        composite = np.zeros(partition_masks[0].shape[1:])
        colors = plt.cm.Set3(np.linspace(0, 1, len(partition_masks)))
        
        for i, mask in enumerate(partition_masks):
            composite[mask[slice_idx]] = i + 1
        
        im = ax.imshow(composite, cmap='Set3', vmin=0, vmax=len(partition_masks))
        ax.set_title(f'Radial Partitions - Slice {slice_idx}')
        ax.axis('off')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Zone ID')
        
    else:
        # 3D visualization
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(partition_masks)))
        

        for i, mask in enumerate(partition_masks):
            coords = np.where(mask)
            if mask.ndim == 3:
                z, y, x = coords
            elif mask.ndim == 2:
                y, x = coords
                z = np.zeros_like(y)  # 用0填充z，实现“平面”可视化
            else:
                raise ValueError(f"Unexpected mask ndim: {mask.ndim}")
            if len(y) > 0:
                # Sample points for visualization (to avoid overcrowding)
                N = len(y)
                if N > 1000:
                    indices = np.random.choice(N, 1000, replace=False)
                    x, y, z = x[indices], y[indices], z[indices]
                ax.scatter(x, y, z, c=[colors[i]], alpha=0.6, s=1, label=f'Zone {i}')


        ax.set_xlabel('X')
        ax.set_ylabel('Y') 
        ax.set_zlabel('Z')
        ax.set_title('3D Radial Partitions')
        ax.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # plt.show()

def plot_partition_features(features_df, save_path=None):
    """
    Plot partition feature analysis results.
    
    Parameters
    ----------
    features_df : pandas.DataFrame
        Feature table from compute_partition_features
    save_path : str, optional
        Path to save the plot
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Zone volume distribution
    axes[0, 0].bar(features_df['zone_id'], features_df['zone_volume'])
    axes[0, 0].set_title('Zone Volume Distribution')
    axes[0, 0].set_xlabel('Zone ID')
    axes[0, 0].set_ylabel('Volume (voxels)')
    
    # Organelle count per zone
    organelle_cols = [col for col in features_df.columns if col.endswith('_count')]
    if organelle_cols:
        for col in organelle_cols:
            org_name = col.replace('_count', '')
            axes[0, 1].plot(features_df['zone_id'].to_numpy(), features_df[col].to_numpy(), 
                           marker='o', label=org_name)


        axes[0, 1].set_title('Organelle Count per Zone')
        axes[0, 1].set_xlabel('Zone ID')
        axes[0, 1].set_ylabel('Count')
        axes[0, 1].legend()
    
    # Organelle density per zone
    density_cols = [col for col in features_df.columns if col.endswith('_density')]
    if density_cols:
        for col in density_cols:
            org_name = col.replace('_density', '')
            axes[1, 0].plot(features_df['zone_id'].to_numpy(), features_df[col].to_numpy(), 
                           marker='s', label=org_name)
        axes[1, 0].set_title('Organelle Density per Zone')
        axes[1, 0].set_xlabel('Zone ID')
        axes[1, 0].set_ylabel('Density (count/volume)')
        axes[1, 0].legend()
    
    # Zone centers distribution
    axes[1, 1].scatter(features_df['zone_center_x'], features_df['zone_center_y'], 
                      c=features_df['zone_id'], cmap='viridis', s=100)
    axes[1, 1].set_title('Zone Centers (X-Y projection)')
    axes[1, 1].set_xlabel('X Center')
    axes[1, 1].set_ylabel('Y Center')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

def visualize_complete_scene(masks, mask_names=None, title="3D Scene", z_scale=1, save=True):
    """
    Simplified 3D scene visualization for 1-4 masks.
    
    Parameters
    ----------
    masks : list of numpy.ndarray
        List of 1-4 mask arrays (e.g., [cell_mask, nucleus_mask, isg_mask, mito_mask])
        All masks should be binary (0/1) data
    mask_names : list of str, optional
        Names for each mask (default: ['Mask 1', 'Mask 2', ...])
    title : str
        Title for the visualization
    z_scale : float
        Z-axis scaling factor (default: 1)
    save : bool
        Whether to save the HTML file (default: True)
    
    Returns
    -------
    plotly.graph_objects.Figure
        The 3D visualization figure
    """
    import plotly.graph_objects as go
    from skimage import measure
    import os
    
    print(f"=== 3D Scene Visualization: {title} ===")
    print(f"Number of masks: {len(masks)}")
    print(f"Z-axis scaling factor: {z_scale}")
    
    # Validate inputs
    if not masks:
        raise ValueError("Please provide masks.")
    
    # Default mask names
    if mask_names is None:
        default_names = ['Cell', 'Nucleus', 'ISG', 'Mitochondria']
        mask_names = default_names[:len(masks)]
    
    # Subsample for performance
    subsample = 3
    masks_ss = []
    for mask in masks:
        if mask is not None:
            masks_ss.append(mask[::subsample, ::subsample, ::subsample])
        else:
            masks_ss.append(None)
    
    print(f"Shape after subsampling: {[m.shape if m is not None else 'None' for m in masks_ss]}")
    
    # Define colors for each mask
    colors = ['red', 'purple', 'blue', 'green']
    opacities = [0.15, 0.25, 0.35, 0.45]
    
    fig = go.Figure()
    
    # Add each mask as a 3D surface
    for i, (mask_ss, name, color, opacity) in enumerate(zip(masks_ss, mask_names, colors, opacities)):
        if mask_ss is None:
            continue
            
        try:
            # Masks are already binary, just ensure float32 for marching cubes
            mask_binary = mask_ss.astype(np.float32)
            
            if np.any(mask_binary > 0):
                vertices, faces, _, _ = measure.marching_cubes(mask_binary, level=0.5)
                
                mesh = go.Mesh3d(
                    x=vertices[:, 2] * subsample,
                    y=vertices[:, 1] * subsample,
                    z=vertices[:, 0] * subsample * z_scale,
                    i=faces[:, 0],
                    j=faces[:, 1],
                    k=faces[:, 2],
                    opacity=opacity,
                    color=color,
                    name=name
                )
                fig.add_trace(mesh)
                print(f"✓ Added {name} surface")
            else:
                print(f"⚠️ {name} mask is empty, skipping")
                
        except Exception as e:
            print(f"✗ Failed to render {name} surface: {e}")
            # Fallback: use scatter plot
            try:
                z_coords, y_coords, x_coords = np.where(mask_binary > 0)
                if len(z_coords) > 2000:
                    indices = np.random.choice(len(z_coords), 2000, replace=False)
                    z_coords, y_coords, x_coords = z_coords[indices], y_coords[indices], x_coords[indices]
                
                scatter = go.Scatter3d(
                    x=x_coords * subsample,
                    y=y_coords * subsample,
                    z=z_coords * subsample * z_scale,
                    mode='markers',
                    marker=dict(size=1, color=color, opacity=opacity*2),
                    name=f'{name} (points)'
                )
                fig.add_trace(scatter)
                print(f"✓ Added {name} as scatter points")
            except Exception as e2:
                print(f"✗ Failed to render {name} as points: {e2}")

    # Update layout
    fig.update_layout(
        title=f"{title} | Z-scale: {z_scale}x",
        width=1200,
        height=900,
        scene=dict(
            xaxis=dict(title="X (pixels)"),
            yaxis=dict(title="Y (pixels)"),
            zaxis=dict(title=f"Z (pixels × {z_scale})" if z_scale != 1 else "Z (pixels)"),
            aspectmode='data',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                up=dict(x=0, y=0, z=1)
            )
        ),
        legend=dict(x=0.02, y=0.98)
    )
    
    # Save and display
    if save:
        # Create results directory if it doesn't exist
        results_dir = os.path.join(os.getcwd(), "results")
        os.makedirs(results_dir, exist_ok=True)
        
        # Generate filename
        safe_title = title.replace(' ', '_').replace('/', '_').replace('\\', '_')
        filename = f"3D_scene_{safe_title}_zscale{z_scale}.html"
        filepath = os.path.join(results_dir, filename)
        
        fig.write_html(filepath)
        print(f"✓ 3D scene saved to: {filepath}")
        
        # Try to open in browser
        try:
            import webbrowser
            webbrowser.open(f'file://{filepath}')
            print("✓ 3D scene opened in browser")
        except:
            print(f"Please open manually in browser: {filepath}")
                
    

    return fig



def plot_radial_rdf(rdf_results, save_path=None):
    """
    Plot radial RDF results
    
    Parameters
    ----------
    rdf_results : dict
        RDF results from calculate_radial_distribution_rdf
    save_path : str, optional
        Path to save the plot
    """
    if rdf_results is None:
        print("No RDF results to plot")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    radii = np.array(rdf_results['radii'])
    rdf_values = np.array(rdf_results['rdf'])
    
    # Plot RDF curve
    ax1.plot(radii, rdf_values, 'b-o', linewidth=2, markersize=6, label='Organelle RDF')
    ax1.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Random distribution')
    ax1.set_xlabel('Normalized Radial Position')
    ax1.set_ylabel('Radial Distribution Function g(r)')
    ax1.set_title('Organelle Radial Distribution Function')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    
    # Plot organelle and cytoplasm counts
    bin_centers = radii
    organelle_counts = rdf_results['organelle_layer_counts']
    cyto_counts = rdf_results['cyto_layer_counts']
    
    width = 0.35
    x = np.arange(len(bin_centers))
    

    ax2.bar(x - width/2, organelle_counts, width, label='Organelle counts', alpha=0.7)
    ax2.bar(x + width/2, cyto_counts, width, label='Cytoplasm counts', alpha=0.7)

    ax2.set_yscale('log')

    ax2.set_xlabel('Radial Bin')
    ax2.set_ylabel('Counts (log scale)')
    ax2.set_title('Organelle and Cytoplasm Distribution')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{r:.2f}' for r in bin_centers])
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"RDF plot saved to: {save_path}")
    
    # plt.show()
    return fig




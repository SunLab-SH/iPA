
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

def calculate_3d_velocity(positions_t1, positions_t2, time_interval_s):
    """
    Calculate 3D velocity from consecutive frame positions.
    
    The 3D velocity is calculated as the Euclidean distance of positions at two 
    consecutive time steps, divided by the time interval between the steps.
    
    Parameters
    ----------
    positions_t1 : np.ndarray
        Organelle positions at time t1 [N, 3] (in micrometers)
    positions_t2 : np.ndarray
        Organelle positions at time t2 [N, 3] (in micrometers)
    time_interval_s : float
        Time interval between frames (seconds)
        
    Returns
    -------
    np.ndarray
        Array of 3D velocities [N] in μm/s
        
    Example
    -------
    >>> import numpy as np
    >>> from ipa.analysis import calculate_3d_velocity
    >>> pos1 = np.array([[10, 20, 30], [15, 25, 35]])  # μm
    >>> pos2 = np.array([[10.5, 20.3, 30.2], [15.4, 25.2, 35.1]])  # μm
    >>> dt = 8.0  # seconds
    >>> velocities = calculate_3d_velocity(pos1, pos2, dt)
    >>> print(velocities)  # [0.064, 0.051] μm/s
    """
    if len(positions_t1) != len(positions_t2):
        raise ValueError("Position arrays must have same length")
    
    if time_interval_s <= 0:
        raise ValueError("Time interval must be positive")
    
    # Calculate Euclidean distances
    displacements = np.linalg.norm(positions_t2 - positions_t1, axis=1)
    
    # Convert to velocity (assuming positions are in μm)
    velocities = displacements / time_interval_s
    
    return velocities


def calculate_radial_velocity(positions, velocities, pm_mask):
    """
    Calculate radial velocity of organelles projected along the direction to the nearest PM voxel.
    
    Parameters
    ----------
    positions : np.ndarray
        Array of organelle positions [N, 3].
    velocities : np.ndarray
        Array of 3D velocity vectors [N, 3].
    pm_mask : np.ndarray
        3D binary mask of the plasma membrane.
        
    Returns
    -------
    np.ndarray
        Array of radial velocity values.
    """
    from scipy.spatial import cKDTree
    
    # Get coordinates of PM voxels
    pm_coords = np.array(np.where(pm_mask > 0)).T
    if len(pm_coords) == 0:
        raise ValueError("PM mask is empty. Cannot calculate radial velocity.")
        
    tree = cKDTree(pm_coords)
    radial_velocities = []
    
    for pos, vel in zip(positions, velocities):
        # Find nearest PM point
        _, idx = tree.query(pos)
        nearest_pm = pm_coords[idx]
        
        # Vector from organelle to PM
        vec_to_pm = nearest_pm - pos
        norm_vec = vec_to_pm / (np.linalg.norm(vec_to_pm) + 1e-8)
        
        # Project velocity onto this radial vector
        v_radial = np.dot(vel, norm_vec)
        radial_velocities.append(v_radial)
        
    return np.array(radial_velocities)


def create_velocity_violin_plot(speed_data, max_speed=30, y_limit=0.35,
                            figure_size=(1.5, 1.2)):
    """
    Create a violin plot for velocity data comparison between categories
    
    Parameters:
    -----------
    speed_data : list
        List of tuples (speed, category)
    max_speed : float, optional
        Maximum speed threshold for filtering outliers
    y_limit : float, optional
        Y-axis upper limit for the plot
    figure_size : tuple, optional
        Figure size as (width, height) in inches
        
    Returns:
    --------
    matplotlib.figure.Figure : The created figure object
    """
    
    # Create DataFrame from speed data
    df = pd.DataFrame(speed_data, columns=['Speed', 'Category'])
    
    # Filter out extreme outliers
    df = df[df['Speed'] <= max_speed]
    
    print(f"Data summary after filtering (speed <= {max_speed}):")
    print(f"Total data points: {len(df)}")
    print(f"LG category: {len(df[df['Category'] == 'LG'])} points")
    print(f"HG category: {len(df[df['Category'] == 'HG'])} points")
    
    # Create figure
    fig = plt.figure(figsize=figure_size, dpi=300)
    
    # Create violin plot
    violin = sns.violinplot(
        x='Category',
        y='Speed',
        data=df,
        palette={'LG': '#00BFBF', 'HG': '#FFA000'},
        inner=None,  # Don't show internal lines
        linewidth=0,  # Remove border lines
        bw=0.2,      # Bandwidth for kernel density estimation
        cut=0,       # Don't extend beyond data range
        scale='area' # Scale violin area
    )
    
    # Add mean lines for each category
    for i, category in enumerate(['LG', 'HG']):
        mean_val = df[df['Category'] == category]['Speed'].mean()
        plt.hlines(mean_val, 
                   xmin=i-0.4, xmax=i+0.4, 
                   colors='black',
                   linewidth=0.8,
                   linestyle='--')
        print(f"Mean velocity for {category}: {mean_val:.4f} um/s")
    
    # Customize plot appearance
    plt.tick_params(axis='x', which='both', length=0)  # Remove x-axis tick marks
    
    # Adjust axes
    plt.ylim(0, y_limit)
    plt.yticks([0.0, 0.2], ['0.0', '0.2'])  # Show only two y-axis ticks
    plt.xlabel('')  
    plt.xticks(ticks=[0, 1], labels=['', ''])  # Remove x-axis labels
    
    # Set axis label styling
    plt.ylabel('3D velocity of ISG [um/s]', labelpad=5)
    plt.tick_params(axis='both', which='major', pad=1)
    
    plt.tight_layout(pad=0.05)
    
    # # Save figure if path is provided
    # if save_path:
    #     plt.savefig(save_path, dpi=900, bbox_inches='tight')
    #     print(f"Figure saved to: {save_path}")
    
    return fig

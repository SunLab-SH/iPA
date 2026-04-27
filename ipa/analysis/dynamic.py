
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns







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

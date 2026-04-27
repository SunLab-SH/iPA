Dynamics Analysis
=================

The dynamics module provides tools for analyzing temporal changes and movement patterns in cellular imaging data.


Functions
---------

Velocity Analysis
~~~~~~~~~~~~~~~~~

.. autofunction:: ipa.analysis.create_velocity_violin_plot

   Creates violin plots for velocity distribution analysis of cellular components.
   
   **Parameters:**
   
   * `velocity_data` (dict or list): Velocity measurements for different conditions or time points
   * `labels` (list): Labels for different datasets
   * `title` (str): Plot title
   * `ylabel` (str): Y-axis label (default: "Velocity")
   * `save_path` (str, optional): Path to save the plot
   
   **Returns:**
   
   * Matplotlib figure object with violin plot

Usage Examples
--------------

Velocity Distribution Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ipa.analysis import create_velocity_violin_plot
   import numpy as np
   
   # Generate example velocity data for different conditions
   velocity_data = {
       'Control': np.random.normal(5.0, 1.5, 100),
       'Treatment_A': np.random.normal(7.2, 2.0, 100), 
       'Treatment_B': np.random.normal(3.8, 1.2, 100)
   }
   
   # Create violin plot
   fig = create_velocity_violin_plot(
       velocity_data,
       labels=['Control', 'Treatment A', 'Treatment B'],
       title='Organelle Velocity Distribution',
       ylabel='Velocity (μm/s)',
       save_path='velocity_analysis.png'
   )

Multi-Condition Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Time-series velocity analysis
   timepoints = ['0h', '2h', '4h', '6h', '8h']
   velocity_by_time = {}
   
   for t in timepoints:
       # Load velocity data for each timepoint
       velocity_by_time[t] = load_velocity_data(f'data_{t}.csv')
   
   # Visualize temporal changes
   fig = create_velocity_violin_plot(
       velocity_by_time,
       labels=timepoints,
       title='Velocity Changes Over Time',
       ylabel='Velocity (μm/s)'
   )


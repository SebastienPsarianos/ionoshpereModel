import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import math


def plot_heatmap_multiple(filename):
    """ 
    Reads coordinate and potential data from a file and plots 2D heatmaps.
    Detects multiple value columns and creates subplots for each.
    """
    # 1. Load the data
    try:
        if not os.path.exists(filename):
            print(f"Error: Could not find the file '{filename}'.")
            return
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Error loading file '{filename}': {e}")
        return

    # 2. Safety Checks
    if data.size == 0:
        print("Error: The data file is empty.")
        return

    # Ensure 2D array even if single line
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # 3. Extract Coordinates and Values
    # Theta is col 0, Phi is col 1.
    # Everything from col 2 onwards are considered values to plot.
    theta_flat = data[:, 0]
    phi_flat = data[:, 1]

    # This selects all columns from index 2 to the end
    values_data = data[:, 2:]

    num_data_points, num_value_cols = values_data.shape

    if num_value_cols == 0:
        print("Error: File contains coordinates but no value columns.")
        return

    print(f"Found {num_value_cols} value column(s) to plot.")

    # 4. Auto-detect Grid Dimensions
    unique_theta = np.unique(theta_flat)
    unique_phi = np.unique(phi_flat)

    n_th = len(unique_theta)
    n_ph = len(unique_phi)

    print(
        f"Detected Grid: nTh={n_th}, nPh={n_ph}(Total rows: {num_data_points})")

    if n_th * n_ph != num_data_points:
        print("Error: Grid dimensions do not match the number of data points.")
        return

    # Create Meshgrid for plotting (converted to degrees)
    # We create this once, as it applies to all subplots
    Phi_grid, Theta_grid = np.meshgrid(
        np.degrees(unique_phi), np.degrees(unique_theta))

    # 5. Setup Dynamic Subplots
    # Determine how many rows/cols for the figure based on num_value_cols
    # We max out at 3 columns wide to keep it readable
    cols_per_row = 3
    n_cols_plot = min(num_value_cols, cols_per_row)
    n_rows_plot = math.ceil(num_value_cols / n_cols_plot)

    fig, axes = plt.subplots(n_rows_plot, n_cols_plot,
                             figsize=(5 * n_cols_plot, 4 * n_rows_plot),
                             constrained_layout=True)

    # If there is only 1 plot, axes is not a list, so we wrap it
    if num_value_cols == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten()

    # 6. Loop through value columns and plot
    for i in range(num_value_cols):
        ax = axes_flat[i]

        # Extract current column and reshape
        current_val_flat = values_data[:, i]
        Z = current_val_flat.reshape((n_th, n_ph))

        # Plot Heatmap
        cp = ax.pcolormesh(Phi_grid, Theta_grid, Z,
                           shading='auto', cmap='viridis')

        # Add colorbar specific to this subplot
        cbar = plt.colorbar(cp, ax=ax)
        cbar.set_label(f'Value (Col {i+2})')

        # Formatting
        ax.set_title(f'Column {i+2}')
        ax.set_xlabel('Phi (deg)')
        ax.set_ylabel('Theta (deg)')
        ax.invert_yaxis()  # North pole at top

    # Hide any unused subplots (if grid is larger than number of plots)
    for j in range(num_value_cols, len(axes_flat)):
        axes_flat[j].axis('off')

    fig.suptitle(
        f'Data Visualization: {os.path.basename(filename)}', fontsize=16)
    plt.show()


# --- Command Line Argument Handling ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot 2D heatmaps for all value columns in a file.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        type=str,
        default='../data/solvedPotential.txt',
        help='Path to the input file.'
    )

    args = parser.parse_args()

    plot_heatmap_multiple(args.input)

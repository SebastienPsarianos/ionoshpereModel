from matplotlib import cm, colors
import matplotlib.pyplot as plt
import math
import os
import argparse
import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')  # Switch to Qt backend


def plot_sphere_multiple(filename):
    """ 
    Reads coordinate and potential data from a file and plots 3D spherical heatmaps.
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
    values_data = data[:, 2:]

    num_data_points, num_value_cols = values_data.shape

    if num_value_cols == 0:
        print("Error: File contains coordinates but no value columns.")
        return

    print(f"Found {num_value_cols} value column(s) to plot.")

    # 4. Auto-detect Grid Dimensions
    # We assume the data forms a regular grid (n_theta x n_phi)
    unique_theta = np.unique(theta_flat)
    unique_phi = np.unique(phi_flat)

    n_th = len(unique_theta)
    n_ph = len(unique_phi)

    print(f"Detected Grid: nTh={n_th}, nPh={
          n_ph} (Total rows: {num_data_points})")

    if n_th * n_ph != num_data_points:
        print("Error: Grid dimensions do not match the number of data points.")
        # Try to proceed anyway if it's just a sorting issue, but warn heavily
        # Often data isn't sorted by theta then phi, so reshaping might scramble it
        # ideally we should re-sort, but let's stick to the previous logic first.

    # 5. Coordinate Conversion (Spherical -> Cartesian)
    # Create the grid mesh
    # Note: meshgrid default is xy indexing (Cartesian), we want ij (Matrix) for easier reshaping
    Theta_grid, Phi_grid = np.meshgrid(unique_theta, unique_phi, indexing='ij')

    # Convert to Cartesian Coordinates for the unit sphere (Radius = 1)
    # Standard Physics Convention (ISO):
    # x = r * sin(theta) * cos(phi)
    # y = r * sin(theta) * sin(phi)
    # z = r * cos(theta)
    R = 1.0
    X = R * np.sin(Theta_grid) * np.cos(Phi_grid)
    Y = R * np.sin(Theta_grid) * np.sin(Phi_grid)
    Z = R * np.cos(Theta_grid)

    # 6. Setup Dynamic Subplots
    cols_per_row = 3
    n_cols_plot = min(num_value_cols, cols_per_row)
    n_rows_plot = math.ceil(num_value_cols / n_cols_plot)

    # Create figure with 3D projection
    fig = plt.figure(figsize=(5 * n_cols_plot, 5 * n_rows_plot),
                     constrained_layout=True)

    # Loop through value columns and plot
    for i in range(num_value_cols):
        # Create subplot with 3d projection
        ax = fig.add_subplot(n_rows_plot, n_cols_plot, i + 1, projection='3d')

        # Extract current column and reshape to match the grid
        current_val_flat = values_data[:, i]

        # NOTE: If your text file is sorted differently (e.g. Phi changes first),
        # this reshape might need to be transposed. Assuming standard C++ nested loop output here.
        V = current_val_flat.reshape((n_th, n_ph))

        # Normalize the data for color mapping (Map values to 0-1)
        norm = colors.Normalize(vmin=V.min(), vmax=V.max())

        # Apply the colormap (viridis) to the normalized values
        # facecolors expects an RGBA array corresponding to the grid
        mapped_colors = cm.viridis(norm(V))

        # Plot the surface
        # rstride/cstride=1 ensures we plot every grid point (high res)
        # shade=False is CRITICAL: We want color to represent Value, not lighting/shadows
        surf = ax.plot_surface(X, Y, Z, facecolors=mapped_colors,
                               rstride=1, cstride=1, shade=False)

        # Create a scalar mappable for the colorbar (since surf doesn't hold the values directly)
        m = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
        m.set_array([])
        cbar = plt.colorbar(m, ax=ax, shrink=0.7, pad=0.1)
        cbar.set_label(f'Value (Col {i+2})')

        # Formatting
        ax.set_title(f'Column {i+2} (Spherical)')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Set equal aspect ratio so the sphere looks round
        ax.set_box_aspect([1, 1, 1])

        # Set initial view angle (Elevation 30, Azimuth 45)
        ax.view_init(elev=30, azim=45)

        # Turn off grid/axis numbers if you want a cleaner look
        # ax.axis('off')

    fig.suptitle(f'3D Visualization: {
                 os.path.basename(filename)}', fontsize=16)
    plt.show()


# --- Command Line Argument Handling ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot 3D spherical heatmaps for all value columns in a file.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        type=str,
        default='../data/solvedPotential.txt',
        help='Path to the input file.'
    )

    args = parser.parse_args()

    plot_sphere_multiple(args.input)

import sys
import numpy as np
import matplotlib.pyplot as plt


def compare_grid_files(file1_path, file2_path):
    # Load data from text files
    try:
        data1 = np.loadtxt(file1_path)
        data2 = np.loadtxt(file2_path)
    except Exception as e:
        print(f"Error reading files: {e}")
        return

    # --- 1. Check Structure ---
    if data1.shape != data2.shape:
        print("File structures do not match.")
        return

    # Extract Theta, Phi, and Values
    theta = data1[:, 0]
    phi = data1[:, 1]
    values1 = data1[:, 2:]
    values2 = data2[:, 2:]

    # Calculate difference
    diff = values1 - values2
    abs_diff = np.abs(diff)

    # --- 2. Calculate Percent Difference ---
    sum_diff = np.sum(abs_diff, axis=0)
    baseline_sum = np.sum(np.abs(values1), axis=0)

    percent_diff = np.zeros_like(sum_diff)
    valid_mask = baseline_sum != 0
    percent_diff[valid_mask] = (
        sum_diff[valid_mask] / baseline_sum[valid_mask]) * 100

    num_value_cols = values1.shape[1]

    # --- 3. Simple Print Output ---
    for i in range(num_value_cols):
        print(f"Value {i+1}: {percent_diff[i]:.4f}%")

    # --- 4. Plotting Heatmaps ---
    # Determine the unique grid values to reshape data into 2D
    unique_theta = np.unique(theta)
    unique_phi = np.unique(phi)
    n_theta = len(unique_theta)
    n_phi = len(unique_phi)

    # Verify if the data forms a perfect grid
    if n_theta * n_phi != len(theta):
        return

    # Reshape coordinates for plotting axes
    Theta, Phi = np.meshgrid(unique_theta, unique_phi, indexing='ij')

    for i in range(num_value_cols):
        # Reshape the 1D difference array into a 2D grid
        grid_diff = diff[:, i].reshape(n_theta, n_phi)

        plt.figure(figsize=(10, 6))

        # Using 'coolwarm' colormap centered at 0
        max_abs_val = np.max(np.abs(grid_diff))
        if max_abs_val == 0:
            max_abs_val = 1e-9  # Prevent warning if diff is perfectly zero

        mesh = plt.pcolormesh(Phi, Theta, grid_diff, cmap='coolwarm', shading='auto',
                              vmin=-max_abs_val, vmax=max_abs_val)

        plt.colorbar(mesh, label='Difference (File 1 - File 2)')
        plt.xlabel('Phi')
        plt.ylabel('Theta')
        plt.title(
            f'Difference Heatmap for Value {i+1} (Total Diff: {percent_diff[i]:.2f}%)')

        filename = f'diff_heatmap_value_{i+1}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare.py <file_1.txt> <file_2.txt>")
    else:
        compare_grid_files(sys.argv[1], sys.argv[2])

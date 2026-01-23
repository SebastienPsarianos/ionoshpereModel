import sys
import numpy as np
import matplotlib.pyplot as plt


def compare_grid_files(file1_path, file2_path):
    print(f"Comparing '{file1_path}' and '{file2_path}'...\n")

    # Load data from text files
    try:
        data1 = np.loadtxt(file1_path)
        data2 = np.loadtxt(file2_path)
    except Exception as e:
        print(f"Error reading files: {e}")
        return

    # --- 1. Check Structure ---
    if data1.shape != data2.shape:
        print("❌ ERROR: File structures do not match.")
        return

    if not np.allclose(data1[:, 0:2], data2[:, 0:2]):
        print("⚠️ WARNING: The theta/phi coordinates do not match exactly.")
    else:
        print("✅ Grid structure and coordinates match perfectly.")

    # Extract Theta, Phi, and Values
    theta = data1[:, 0]
    phi = data1[:, 1]
    values1 = data1[:, 2:]
    values2 = data2[:, 2:]

    # Calculate difference
    diff = values1 - values2
    abs_diff = np.abs(diff)

    # --- 2 & 3. Calculate Sum and Normalize ---
    sum_diff = np.sum(abs_diff, axis=0)
    baseline_sum = np.sum(np.abs(values1), axis=0)

    percent_diff = np.zeros_like(sum_diff)
    valid_mask = baseline_sum != 0
    percent_diff[valid_mask] = (
        sum_diff[valid_mask] / baseline_sum[valid_mask]) * 100

    num_value_cols = values1.shape[1]

    print("\n" + "-" * 55)
    print(f"{'Value Col':<12} | {'Sum of Abs Diff':<18} | {
          'Percent Diff (%)':<15}")
    print("-" * 55)
    for i in range(num_value_cols):
        print(f"Value {i+1:<6} | {sum_diff[i]
              :<18.6f} | {percent_diff[i]:<15.4f}%")
    print("-" * 55)

    # --- 4. Plotting Heatmaps ---
    print("\nGenerating heatmaps...")

    # Determine the unique grid values to reshape data into 2D
    unique_theta = np.unique(theta)
    unique_phi = np.unique(phi)

    n_theta = len(unique_theta)
    n_phi = len(unique_phi)

    # Verify if the data forms a perfect grid
    if n_theta * n_phi != len(theta):
        print("⚠️ Cannot generate heatmaps: Data is not on a regular rectangular grid.")
        return

    # Reshape coordinates for plotting axes
    Theta, Phi = np.meshgrid(unique_theta, unique_phi, indexing='ij')

    for i in range(num_value_cols):
        # Reshape the 1D difference array into a 2D grid
        grid_diff = diff[:, i].reshape(n_theta, n_phi)

        # Create the plot
        plt.figure(figsize=(10, 6))

        # Using 'coolwarm' colormap to show positive/negative differences clearly
        # Centering the colormap at 0 using vmin and vmax
        max_abs_val = np.max(np.abs(grid_diff))
        mesh = plt.pcolormesh(Phi, Theta, grid_diff, cmap='coolwarm', shading='auto',
                              vmin=-max_abs_val, vmax=max_abs_val)

        plt.colorbar(mesh, label='Difference (File 1 - File 2)')
        plt.xlabel('Phi')
        plt.ylabel('Theta')
        plt.title(f'Difference Heatmap for Value {
                  i+1} (Total Diff: {percent_diff[i]:.2f}%)')

        # Save the plot
        filename = f'diff_heatmap_value_{i+1}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Saved plot: {filename}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_files.py <file_1.txt> <file_2.txt>")
    else:
        compare_grid_files(sys.argv[1], sys.argv[2])

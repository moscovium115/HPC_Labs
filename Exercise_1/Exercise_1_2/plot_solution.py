import numpy as np
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_and_plot_poisson():
    # 1. Find all output files
    file_pattern = os.path.join("data", "output*.dat")
    files = glob.glob(file_pattern)

    if not files:
        print("No output files found in 'data/'. Did you run the MPI code?")
        return

    print(f"Found {len(files)} output files. Loading data...")

    # 2. Load data from all files
    all_data = []
    for fname in files:
        try:
            data_chunk = np.loadtxt(fname)
            
            if data_chunk.size == 0:
                continue
            if data_chunk.ndim == 1:
                data_chunk = data_chunk.reshape(1, -1)
                
            all_data.append(data_chunk)
        except Exception as e:
            print(f"Warning: Could not read {fname}: {e}")

    if not all_data:
        print("No valid data found.")
        return

    # Combine all chunks
    combined_data = np.vstack(all_data)

    # 3. Determine grid dimensions
    # CHANGE 1: We treat the C-code X (col 0) as horizontal width
    #           and C-code Y (col 1) as vertical height.
    max_x_idx = int(combined_data[:, 0].max()) 
    max_y_idx = int(combined_data[:, 1].max())

    # CHANGE 2: Grid shape is now (Height, Width) -> (Y, X)
    grid = np.zeros((max_y_idx + 2, max_x_idx + 2))

    # 4. Fill the grid
    # CHANGE 3: Swap indices. grid[row, col] becomes grid[y, x]
    # We use the raw file columns: col 0 is X, col 1 is Y
    raw_x = combined_data[:, 0].astype(int)
    raw_y = combined_data[:, 1].astype(int)
    values = combined_data[:, 2]

    # Assign values: grid[y_index, x_index] = value
    grid[raw_y, raw_x] = values

    # 5. Plotting
    plt.figure(figsize=(10, 8))
    
    # CHANGE 4: 
    # - origin='lower': puts (0,0) at bottom-left (Standard Cartesian)
    # - updated labels to match your request
    plt.imshow(grid, cmap='inferno', origin='lower', interpolation='nearest')
    
    plt.colorbar(label='Phi (Potential)')
    plt.title(f'2D Poisson Solution\nGrid Size: {max_x_idx} x {max_y_idx}')
    
    # X-axis is now the horizontal axis
    plt.xlabel('X Index (Horizontal)')
    # Y-axis is now the vertical axis
    plt.ylabel('Y Index (Vertical)')
    
    plt.savefig('poisson_solution.png', dpi=150)
    print("Plot saved as 'poisson_solution.png'")

if __name__ == "__main__":
    load_and_plot_poisson()
import re
import pandas as pd
import numpy as np

# content of the file to parse (reads from local file 'ex4_2_output.txt')
filename = 'ex4_2_output.txt'

def parse_mpi_output(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # Split content by the separator lines to handle different grid sizes
    blocks = content.split('----------------------------------------')
    
    data_list = []

    # Iterate through blocks to find grid sections
    current_grid = None
    
    for block in blocks:
        # Check for grid size header
        grid_match = re.search(r'Running simulation for Grid:\s+(\d+\s+x\s+\d+)', block)
        if grid_match:
            current_grid = grid_match.group(1).replace(' ', '')
            continue # Move to the data part of this block (or next split)

        if not current_grid:
            continue

        # Check if this block contains the timing data
        if "Time Breakdown" in block:
            # Extract Iterations
            iter_match = re.search(r'Number of iterations\s+:\s+(\d+)', block)
            iterations = int(iter_match.group(1)) if iter_match else 0

            # Extract all timing breakdown lines
            # Updated Format: ... | Idle: 0.0050s
            pattern = r'Calc:\s+([\d\.]+)s\s+\|\s+Exchange:\s+([\d\.]+)s\s+\|\s+Global:\s+([\d\.]+)s\s+\|\s+Idle:\s+([\d\.]+)s'
            matches = re.findall(pattern, block)
            
            if matches:
                # Convert strings to floats
                timings = np.array(matches, dtype=float)
                
                # Calculate means across all processes (ranks)
                avg_calc = np.mean(timings[:, 0])
                avg_exch = np.mean(timings[:, 1])
                avg_glob = np.mean(timings[:, 2])
                avg_idle = np.mean(timings[:, 3]) # New column
                
                # Calculate Ratio (Computation / Communication)
                # Note: Communication cost is Exchange + Global (Idle is usually excluded from pure Comm cost in this ratio)
                total_comm = avg_exch + avg_glob
                ratio = avg_calc / total_comm if total_comm > 0 else 0

                data_list.append({
                    'Grid Size': current_grid,
                    'Iterations': iterations,
                    'Computation (s)': round(avg_calc, 4),
                    'Exchange (s)': round(avg_exch, 4),
                    'Global Comm (s)': round(avg_glob, 4),
                    'Idle (s)': round(avg_idle, 4), # Added to DataFrame
                    'Ratio (Comp/Comm)': round(ratio, 1)
                })
            
            # Reset current grid after processing data to avoid duplicates
            current_grid = None

    return pd.DataFrame(data_list)

# --- Execution ---
try:
    df = parse_mpi_output(filename)
    
    if not df.empty:
        print("\n=== Generated DataFrame ===")
        print(df)
        
        print("\n=== LaTeX Code ===")
        print(df.to_latex(index=False, caption="Performance Breakdown by Grid Size ($P=4$)", label="tab:perf_metrics"))
    else:
        print("No data found. Make sure ex4_2_output.txt exists and is populated.")

except FileNotFoundError:
    print(f"Error: Could not find file '{filename}'. Make sure it is in the same directory.")
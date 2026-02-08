import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# 1. Read Data
df = pd.read_csv('data_ex_1_2_11/comm_analysis.csv')

# 2. Calculate Message Size per Exchange
# Note: The time in CSV is the TOTAL accumulated time for 1000 iterations (2000 calls).
# We need the time per SINGLE call.
df['TimePerCall_us'] = (df['AvgCommTime'] / 2000.0) * 1e6  # Convert seconds to microseconds

# Define message size (L) in Bytes based on topology
def get_message_size(row):
    if row['Topology'] == '4x1':
        return row['GridSize'] * 8
    elif row['Topology'] == '2x2':
        return (row['GridSize'] / 2) * 8
    return 0

df['MessageSize_B'] = df.apply(get_message_size, axis=1)

# 3. Perform Linear Regression for 4x1 (Primary Interest)
subset_4x1 = df[df['Topology'] == '4x1']

# Check if subset is not empty
if len(subset_4x1) > 1:
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        subset_4x1['MessageSize_B'], 
        subset_4x1['TimePerCall_us']
    )
    latency = intercept
    bandwidth_MBps = (1 / slope) if slope != 0 else 0  # Avoid division by zero
    
    print("-" * 40)
    print(f"Analysis for 4x1 Topology (Contiguous Data)")
    print("-" * 40)
    print(f"Latency (alpha):   {latency:.4f} us")
    print(f"Bandwidth (1/beta): {bandwidth_MBps:.2f} MB/s")
    print(f"R-squared:         {r_value**2:.4f}")
    print("-" * 40)
else:
    print("Insufficient data for 4x1 topology")
    slope, intercept = None, None

# 4. Perform Linear Regression for 2x2 (Secondary Interest)
subset_2x2 = df[df['Topology'] == '2x2']

if len(subset_2x2) > 1:
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(
        subset_2x2['MessageSize_B'], 
        subset_2x2['TimePerCall_us']
    )
    print(f"Analysis for 2x2 Topology (Strided Data)")
    print("-" * 40)
    print(f"Latency (alpha):   {intercept2:.4f} us")
    print(f"Bandwidth (1/beta): {(1/slope2):.2f} MB/s" if slope2 != 0 else "Bandwidth: undefined")
    print(f"R-squared:         {r_value2**2:.4f}")
    print("-" * 40)
else:
    print("Insufficient data for 2x2 topology")
    slope2, intercept2 = None, None

# 5. Plotting
plt.figure(figsize=(10, 6))

# Plot 4x1
if len(subset_4x1) > 0:
    plt.plot(subset_4x1['MessageSize_B'], subset_4x1['TimePerCall_us'], 
             'o-', label='4x1 (Contiguous)', markersize=6)
    
    # Plot Regression Line 4x1
    if slope is not None and intercept is not None:
        x_min = subset_4x1['MessageSize_B'].min()
        x_max = subset_4x1['MessageSize_B'].max()
        x_vals = np.linspace(max(0, x_min - 1000), x_max + 1000, 100)
        y_vals = intercept + slope * x_vals
        bandwidth_MBps = (1 / slope) if slope != 0 else 0
        plt.plot(x_vals, y_vals, '--', color='blue', alpha=0.5, 
                label=f'Fit 4x1: {intercept:.2f}us + L/{bandwidth_MBps:.0f} MB/s')

# Plot 2x2
if len(subset_2x2) > 0:
    plt.plot(subset_2x2['MessageSize_B'], subset_2x2['TimePerCall_us'], 
             's-', label='2x2 (Strided)', markersize=6)
    
    # Plot Regression Line 2x2
    if slope2 is not None and intercept2 is not None:
        x_min2 = subset_2x2['MessageSize_B'].min()
        x_max2 = subset_2x2['MessageSize_B'].max()
        x_vals2 = np.linspace(max(0, x_min2 - 1000), x_max2 + 1000, 100)
        y_vals2 = intercept2 + slope2 * x_vals2
        bandwidth2_MBps = (1 / slope2) if slope2 != 0 else 0
        plt.plot(x_vals2, y_vals2, '--', color='orange', alpha=0.5, 
                label=f'Fit 2x2: {intercept2:.2f}us + L/{bandwidth2_MBps:.0f} MB/s')

plt.title('Communication Time vs Message Size')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Time per Call (us)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig('comm_analysis_plot.png', dpi=150)
print("\nPlot saved to comm_analysis_plot.png")
plt.show()
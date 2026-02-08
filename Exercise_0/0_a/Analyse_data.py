import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def parse_mpi_file_bytes(file_path):
    elements, times = [], []
    try:
        with open(file_path, 'r') as f:
            current_n = None
            for line in f:
                send_match = re.search(r'Sending (\d+) elements', line)
                if send_match:
                    current_n = int(send_match.group(1))
                time_match = re.search(r'Ping Pong took ([\d\.e-]+) seconds', line)
                if time_match and current_n is not None:
                    times.append(float(time_match.group(1)))
                    elements.append(current_n)
                    current_n = None
    except FileNotFoundError:
        print(f"Warning: {file_path} not found.")
        return None
    
    els, ts = np.array(elements), np.array(times)
    size_bytes = els * 4
    
    df = pd.DataFrame({
        'size_bytes': size_bytes,
        'communication_time': ts 
    })
    
    return df

# Loading in all data files
df_std_same = parse_mpi_file_bytes('0_a_data/ping_std_same.txt')
df_std_diff = parse_mpi_file_bytes('0_a_data/ping_std_diff.txt')
df_extra_same = parse_mpi_file_bytes('0_a_data/ping_same_extra.txt')
df_extra_diff = parse_mpi_file_bytes('0_a_data/ping_diff_extra.txt')


fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# First subplot: Standard
axes[0].scatter(df_std_same["size_bytes"], df_std_same["communication_time"], 
                alpha=0.6, s=30, label="1 Node (Shared)")
axes[0].scatter(df_std_diff["size_bytes"], df_std_diff["communication_time"], 
                alpha=0.6, s=30, label="2 Nodes (Network)")
axes[0].set_xscale('log', base=2)
axes[0].set_yscale('log')

powers = np.arange(1, 23)
tick_values = 2**powers
tick_labels = [f'$2^{{{p}}}$' for p in powers]
axes[0].set_xticks(tick_values)
axes[0].set_xticklabels(tick_labels, rotation=45)
axes[0].set_xlim(2**0, 2**23)
axes[0].set_xlabel('Message Size (Bytes)')
axes[0].set_ylabel('Communication Time $t_{comm}$ (s)')
axes[0].set_title('MPI_Send/MPI_Recv Performance')
axes[0].legend()
axes[0].grid(True, which='both', linestyle=':', alpha=0.6)

# Second subplot: Extra
axes[1].scatter(df_extra_same["size_bytes"], df_extra_same["communication_time"], 
                alpha=0.6, s=30, label="1 Node (Shared)")
axes[1].scatter(df_extra_diff["size_bytes"], df_extra_diff["communication_time"], 
                alpha=0.6, s=30, label="2 Nodes (Network)")
axes[1].set_xscale('log', base=2)
axes[1].set_yscale('log')
axes[1].set_xticks(tick_values)
axes[1].set_xticklabels(tick_labels, rotation=45)
axes[1].set_xlim(2**0, 2**23)
axes[1].set_xlabel('Message Size (Bytes)')
axes[1].set_ylabel('Communication Time $t_{comm}$ (s)')
axes[1].set_title('MPI_SendRecv MPI Performance')
axes[1].legend()
axes[1].grid(True, which='both', linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('mpi_performance_loglog.png')
plt.show()

#############################

def calculate_bandwidth(df):
    """
    Calculates Bandwidth in MB/s.
    Bandwidth = (Size in Bytes) / (Time in Seconds) / 1,000,000
    """
    df = df.copy()
    # Avoid division by zero if time is 0 (unlikely but safe)
    df['bandwidth_MBs'] = (df['size_bytes'] / df['communication_time']) / 1e6
    return df

# Calculate bandwidth for all datasets
df_std_same = calculate_bandwidth(df_std_same)
df_std_diff = calculate_bandwidth(df_std_diff)
df_extra_same = calculate_bandwidth(df_extra_same)
df_extra_diff = calculate_bandwidth(df_extra_diff)

# Create the plots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# --- Plot 1: Bandwidth for Standard Send/Recv ---
axes[0].plot(df_std_same["size_bytes"], df_std_same["bandwidth_MBs"], 
             marker='o', linestyle='-', label="1 Node (Shared)", alpha=0.8)
axes[0].plot(df_std_diff["size_bytes"], df_std_diff["bandwidth_MBs"], 
             marker='s', linestyle='--', label="2 Nodes (Network)", alpha=0.8)

axes[0].set_xscale('log', base=2)
# We usually use linear scale for Y-axis in bandwidth to see saturation clearly, 
# but log can help if the gap is huge. Let's try Linear first for clarity.
axes[0].set_yscale('linear') 

axes[0].set_title('Bandwidth - Standard Send/Recv')
axes[0].set_ylabel('Bandwidth (MB/s)')
axes[0].set_xlabel('Message Size (Bytes)')
axes[0].grid(True, which="both", linestyle=':', alpha=0.6)
axes[0].legend()

# --- Plot 2: Bandwidth for SendRecv ---
axes[1].plot(df_extra_same["size_bytes"], df_extra_same["bandwidth_MBs"], 
             marker='o', linestyle='-', label="1 Node (Shared)", alpha=0.8)
axes[1].plot(df_extra_diff["size_bytes"], df_extra_diff["bandwidth_MBs"], 
             marker='s', linestyle='--', label="2 Nodes (Network)", alpha=0.8)

axes[1].set_xscale('log', base=2)
axes[1].set_yscale('linear')
axes[1].set_title('Bandwidth - MPI_SendRecv')
axes[1].set_ylabel('Bandwidth (MB/s)')
axes[1].set_xlabel('Message Size (Bytes)')
axes[1].grid(True, which="both", linestyle=':', alpha=0.6)
axes[1].legend()

# --- Formatting Ticks ---
powers = np.arange(1, 23)
tick_values = 2**powers
tick_labels = [f'$2^{{{p}}}$' for p in powers]

for ax in axes:
    ax.set_xticks(tick_values)
    ax.set_xticklabels(tick_labels, rotation=45)
    ax.set_xlim(2**0, 2**23)

plt.tight_layout()
plt.savefig('mpi_bandwidth_analysis.png')
plt.show()

########

# now remove the first row of every dataframe
df_std_same = df_std_same.iloc[1:].reset_index(drop=True)
df_std_diff = df_std_diff.iloc[1:].reset_index(drop=True)
df_extra_same = df_extra_same.iloc[1:].reset_index(drop=True)
df_extra_diff = df_extra_diff.iloc[1:].reset_index(drop=True)
print("First row removed from all dataframes. Dataframes are now:")
print("Standard Send/Recv - Same Node:\n", df_std_same.head())


###########################################


def fit_and_plot_regimes(ax, df, df_name, color, regions):
    """
    Fits Weighted Linear Regressions (minimizing relative error) to specific regions.
    """
    print(df_name)

    for (start_bytes, end_bytes) in regions:
        # Filter data
        subset = df[(df['size_bytes'] >= start_bytes) & (df['size_bytes'] <= end_bytes)]
        
        if len(subset) < 2:
            continue

        X = subset['size_bytes'].values.reshape(-1, 1)
        y = subset['communication_time'].values
        
        # --- NEW: Weighted Least Squares ---
        # Weight by 1/y^2. This balances the influence of small vs large points.
        # It minimizes sum((y_pred - y_true)/y_true)^2 -> Relative Error
        weights = 1.0 / (y**2)
        
        model = LinearRegression()
        # Pass weights here
        model.fit(X, y, sample_weight=weights)
        # -----------------------------------
        
        # Create prediction line (handle log-space plotting start)
        plot_start = start_bytes
        if plot_start <= 0:
            plot_start = max(1, subset['size_bytes'].min()) 
            
        X_pred = np.geomspace(plot_start, end_bytes, 100).reshape(-1, 1)
        y_pred = model.predict(X_pred)
        
        # Plot the fit line
        ax.plot(X_pred, y_pred, color='black', linestyle='--', linewidth=1.5, alpha=0.9)
        
        # Calculate L and B for label
        intercept = model.intercept_
        slope = model.coef_[0]
        bandwidth = 1.0/slope if slope > 0 else 0
        
        # Print slope to console as requested
        print(f"Region {start_bytes}-{end_bytes} | Slope: {slope} | Intercept: {intercept}| Bandwidth: {bandwidth/1e6:.2f} MB/s | R^2: {model.score(X, y):.4f}")

        # Place label near the geometric midpoint
        mid_point = np.sqrt(plot_start * end_bytes)
        y_mid = model.predict([[mid_point]])[0]
        
        # Adjust label position slightly based on regime
        label_y_pos = y_mid * 1.3
        
        # ax.text(mid_point, label_y_pos, 
        #         f'L={intercept*1e6:.1f}µs\nB={bandwidth/1e6:.0f}MB/s', 
        #         fontsize=8, color=color, ha='center', va='bottom',
        #         bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
# --- 3. Plotting Setup ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Define Ticks
powers = np.arange(1, 23)
tick_values = 2.0**powers
tick_labels = [f'$2^{{{int(p)}}}$' for p in powers]

# ==========================================
# LEFT PLOT: Standard Send/Recv
# ==========================================
ax = axes[0]
# Scatter Raw Data
ax.scatter(df_std_same["size_bytes"], df_std_same["communication_time"], label="1 Node", alpha=0.5)
ax.scatter(df_std_diff["size_bytes"], df_std_diff["communication_time"], label="2 Nodes", alpha=0.8)

# Define Regions for Standard (Jump at 2^18)
# Region 1: Eager (0 to 2^17)
# Region 2: Rendezvous (2^18 to Max)
regions_std_1_node = [
    [2**3,2**6],     # Eager
    # [2**7, 2**22],  # Rendezvous
    [2**8, 2**17],
    [2**18, 2**22]

]
# checkpoint

regions_std_2_nodes = [
    [2**3, 2**7],     # Eager
    [2**9, 2**17],  # Rendezvous
    [2**18, 2**22]  # Rendezvous
]

# Fit and Plot lines for Same Node (1 Node) data
fit_and_plot_regimes(ax, df_std_same, "std_same", 'tab:blue', regions_std_1_node)

# Fit and Plot lines for Network (2 Nodes) data
fit_and_plot_regimes(ax, df_std_diff, "std_diff",'tab:orange', regions_std_2_nodes)

ax.set_title('Standard Send/Recv (Log-Log Fit)')
ax.set_xlabel('Message Size (Bytes)')
ax.set_ylabel('Time (s)')
ax.set_xscale('log', base=2)
ax.set_yscale('log')
ax.set_xticks(tick_values)
ax.set_xticklabels(tick_labels, rotation=45)
ax.set_xlim(2**0, 2**23)
ax.legend()
ax.grid(True, which='both', linestyle=':', alpha=0.4)


# ==========================================
# RIGHT PLOT: MPI_SendRecv
# ==========================================
ax = axes[1]
# Scatter Raw Data
ax.scatter(df_extra_same["size_bytes"], df_extra_same["communication_time"], label="1 Node", alpha=0.5)
ax.scatter(df_extra_diff["size_bytes"], df_extra_diff["communication_time"], label="2 Nodes", alpha=0.8)

# Define Regions for SendRecv (Jump earlier at 2^17)
# Region 1: Eager (0 to 2^16)
# Region 2: Rendezvous (2^17 to Max)
regions_extra_node_1 = [
    # [2**3, 2**3],     # Eager
    [2**3,2**6],
    [2**8, 2**17],  # Rendezvous

    [2**18, 2**22]
    # [2**18, 2**22]  # Rendezvous
]

regions_extra_nodes_2 = [
    [2**3, 2**7],     # Eager
    # [2**8, 2**18],  # Rendezvous
    [2**9, 2**17],  # Eager
    [2**18, 2**22]  # Rendezvous
]

# print(df_extra_same.tail(10))
# Fit and Plot lines for Same Node (1 Node) data
fit_and_plot_regimes(ax, df_extra_same, "extra_same",'tab:blue', regions_extra_node_1)

# Fit and Plot lines for Network (2 Nodes) data
fit_and_plot_regimes(ax, df_extra_diff, "extra_diff", 'tab:orange', regions_extra_nodes_2)

ax.set_title('MPI_SendRecv (Log-Log Fit)')
ax.set_xlabel('Message Size (Bytes)')
ax.set_ylabel('Time (s)')
ax.set_xscale('log', base=2)
ax.set_yscale('log')
ax.set_xticks(tick_values)
ax.set_xticklabels(tick_labels, rotation=45)
ax.set_xlim(2**0, 2**23)
ax.legend()
ax.grid(True, which='both', linestyle=':', alpha=0.4)

plt.tight_layout()
plt.savefig('mpi_loglog_fits.png')
plt.show()


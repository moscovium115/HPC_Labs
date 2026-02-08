import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# --- CONFIGURATION ---
INPUT_FILE = 'data_ex_1_2_3/scaling_trace.csv'
OUTPUT_PLOT_SCALING = 'plot_linear_scaling.png'
OUTPUT_PLOT_BETA = 'plot_beta_comparison.png'

def analyze_data():
    # 1. Load Data
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found. Run the bash script first.")
        return

    df = pd.read_csv(INPUT_FILE)
    
    # Sort for consistent plotting
    df = df.sort_values(by=['Grid', 'Topology', 'Iterations'])

    # 2. Setup Regression Storage
    reg_results = []

    # 3. Setup Plotting (Linear Scaling)
    grids = df['Grid'].unique()
    grids.sort()
    
    # Create subplots: One for each grid size
    fig, axes = plt.subplots(1, len(grids), figsize=(6 * len(grids), 5), sharey=False)
    if len(grids) == 1: axes = [axes] # Handle single grid case

    sns.set_style("whitegrid")
    colors = sns.color_palette("bright")

    print(f"{'Grid':<10} {'Topology':<10} {'Alpha (s)':<15} {'Beta (s/iter)':<15} {'R^2':<10}")
    print("-" * 65)

    # 4. Loop through Grids
    for i, grid in enumerate(grids):
        ax = axes[i]
        subset = df[df['Grid'] == grid]
        
        # Plot Scatter Points
        sns.scatterplot(
            data=subset, x='Iterations', y='Time', 
            hue='Topology', style='Topology', 
            s=80, alpha=0.7, ax=ax, palette='bright'
        )

        # Calculate and Plot Regression Lines
        unique_topos = subset['Topology'].unique()
        for j, topo in enumerate(unique_topos):
            topo_data = subset[subset['Topology'] == topo]
            
            # Linear Regression: t = alpha + beta * n
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                topo_data['Iterations'], topo_data['Time']
            )
            
            # Save results
            reg_results.append({
                'Grid': grid, 
                'Topology': topo, 
                'Alpha': intercept, 
                'Beta': slope, 
                'R2': r_value**2
            })
            
            # Print to console
            print(f"{grid:<10} {topo:<10} {intercept:<15.6f} {slope:<15.6f} {r_value**2:<10.4f}")

            # Plot Line
            x_vals = topo_data['Iterations']
            y_vals = intercept + slope * x_vals
            # Find the color corresponding to this topology
            color = colors[j % len(colors)] 
            ax.plot(x_vals, y_vals, '--', linewidth=1, color='black', alpha=0.5)

        ax.set_title(f"Grid: {grid}x{grid}")
        ax.set_ylabel("Time (s)")
        ax.set_xlabel("Iterations")
        ax.legend(title='Topology')

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_SCALING)
    print(f"\nScaling plot saved to {OUTPUT_PLOT_SCALING}")

    # 5. Plot Beta Comparison (Bar Chart)
    # This clearly shows which topology is faster per iteration
    results_df = pd.DataFrame(reg_results)
    
    plt.figure(figsize=(10, 6))
    barplot = sns.barplot(
        data=results_df, 
        x='Grid', 
        y='Beta', 
        hue='Topology',
        palette='viridis'
    )
    
    plt.title("Comparison of Beta (Time per Iteration)")
    plt.ylabel("Beta (seconds/iteration)")
    plt.xlabel("Grid Size")
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_BETA)
    print(f"Beta comparison plot saved to {OUTPUT_PLOT_BETA}")

if __name__ == "__main__":
    analyze_data()
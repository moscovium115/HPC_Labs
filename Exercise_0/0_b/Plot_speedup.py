import pandas as pd
import matplotlib.pyplot as plt
import os

csv_file = "mm_results_final.csv"
if not os.path.exists(csv_file):
    print(f"Error: {csv_file} not found. Run Process_results_v3.py first.")
    exit(1)

# Load data
df = pd.read_csv(csv_file)
plt.figure(figsize=(10, 6))

# Plot Ideal Speedup (Linear Reference)

plt.plot(df["P"], df["P"], 'k--', label="Ideal Linear Speedup ($S = P$)", alpha=0.5)

# Plot Actual Speedup Line
plt.plot(df["P"], df["Speedup"], 'b-', zorder=1, label="Actual Speedup")

# Filter data
shared = df[df["P"] <= 40]
distrib = df[df["P"] > 40]

# Plot Shared Memory (Green dots)
plt.scatter(shared["P"], shared["Speedup"], color='green', s=100, zorder=2, 
            label="Shared Memory (1 Node)")

# Plot Distributed Memory (Red dots)
plt.scatter(distrib["P"], distrib["Speedup"], color='red', s=100, zorder=2, 
            label="Distributed Memory (2 Nodes)")


# 5. LaTeX Labels
plt.xlabel(r"Number of Processes ($P$)", fontsize=12)
plt.ylabel(r"Speedup ($S = T_{serial} / T_{parallel}$)", fontsize=12)
plt.title(r"Matrix Multiplication Scaling: Shared vs. Distributed Memory", fontsize=14)
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Save plot
plt.savefig("speedup_plot.png", dpi=300)
print("Plot saved to speedup_plot.png")
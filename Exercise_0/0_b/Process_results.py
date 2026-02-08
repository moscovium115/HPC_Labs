import pandas as pd
import re
import os

filename = "0_b_data/mm_out_all.txt"
if not os.path.exists(filename):
    filename = "mm_out_all.txt" 

data = []

# Regex to capture P and Time
pattern = re.compile(r"P=(\d+),\s+Time=([0-9\.]+)\s+s")

# Hardcoded layouts based on your SLURM logs
# We use generic "NodeA/NodeB" to keep the table clean
layout_map = {
    41: "40(NodeA) + 1(NodeB)",
    48: "40(NodeA) + 8(NodeB)",
    64: "40(NodeA) + 24(NodeB)"
}

with open(filename, 'r') as f:
    for line in f:
        match = pattern.search(line)
        if match:
            p_val = int(match.group(1))
            time_val = float(match.group(2))
            
            # Determine Layout String
            if p_val <= 40:
                layout_desc = f"{p_val}(NodeA)"
                mode = "Shared Mem"
            else:
                # Use our manual map, or default to generic distributed if unknown
                layout_desc = layout_map.get(p_val, "Distributed")
                mode = "Distributed"
            
            data.append({
                "P": p_val,
                "Time": time_val,
                "Mode": mode,
                "Task_Distribution": layout_desc
            })

# Create DataFrame
df = pd.DataFrame(data)
df.set_index("P", inplace=True)

# Metrics
t_serial = df.loc[1, "Time"]
df["Speedup"] = t_serial / df["Time"]
df["Efficiency"] = df["Speedup"] / df.index

# Reorder columns
df = df[["Time", "Speedup", "Efficiency", "Mode", "Task_Distribution"]]

print("-" * 80)
print("Final Results Table")
print("-" * 80)
print(df)

df.to_csv("mm_results_final.csv")
print("\nSaved to 'mm_results_final.csv'")

latex_table = df.to_latex(
    index=False, 
    caption="Performance Metrics (Sorted by P)", 
    label="tab:performance_results",
    column_format="rcccc l",
    escape=False
)

print(latex_table)
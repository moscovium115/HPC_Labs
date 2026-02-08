import pandas as pd

# Load the data from the csv file
df = pd.read_csv('mm_results_final.csv')

# Ensure the data is sorted by P (increasing order)
df = df.sort_values(by='P', ascending=True)

# Optional: Format numerical columns for better LaTeX presentation
# Using 6 decimal places for Time and 4 for Speedup/Efficiency
df['Time'] = df['Time'].map(lambda x: f'{x:.6f}')
df['Speedup'] = df['Speedup'].map(lambda x: f'{x:.4f}')
df['Efficiency'] = df['Efficiency'].map(lambda x: f'{x:.4f}')

# Generate the LaTeX table code
# We use column_format to align P, Mode, and Distribution appropriately
latex_table = df.to_latex(
    index=False, 
    caption="Performance Metrics (Sorted by P)", 
    label="tab:performance_results",
    column_format="rcccc l",
    escape=False
)

print(latex_table)
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DATA_FILE = 'data_ex_1_2_6/error_trace_800.csv'
OUTPUT_IMG = 'plot_error_convergence.png'

def plot_error():
    try:
        df = pd.read_csv(DATA_FILE)
    except FileNotFoundError:
        print("Data file not found!")
        return

    plt.figure(figsize=(10, 6))
    
    # Log-Linear Plot
    plt.semilogy(df['Iteration'], df['Error'], 'b-', linewidth=2, label='Global Max Error')
    
    # Add a reference line for the tolerance
    plt.axhline(y=1e-4, color='r', linestyle='--', label='Tolerance (1e-4)')

    plt.title('Convergence Behavior (800x800, Omega=1.95)')
    plt.xlabel('Iteration')
    plt.ylabel('Global Maximum Error (log scale)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(OUTPUT_IMG)
    print(f"Plot saved to {OUTPUT_IMG}")

if __name__ == "__main__":
    plot_error()
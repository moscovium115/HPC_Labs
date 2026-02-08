import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# PART 1: PARSING LOGIC
# ==========================================
def parse_lab_results(filename):
    """
    Parses lab_results.out to extract performance data.
    """
    data = []
    
    # Regex to extract the specific timing lines
    rx_cpu = re.compile(r"CPU: run time =\s+([0-9\.]+) secs")
    
    # Robust Regex for GPU times
    # Captures "Unified", "Shared", or "Global" and their respective times
    rx_gpu = re.compile(
        r"GPU (?P<type>Unified|Shared|Global) \(Compute Only\):\s+(?P<compute>[0-9\.]+).*?"
        r"GPU (?P=type) \(Total \+ Mem\):\s+(?P<total>[0-9\.]+)", 
        re.DOTALL
    )

    try:
        with open(filename, 'r') as f:
            content = f.read()
            
        # Split content by the Header to isolate configurations
        segments = content.split(">>> RUNNING:")[1:] 
        
        for segment in segments:
            # 1. Get the Matrix/Block config
            header_match = re.match(r"\s*Matrix Size (\d+) \| Block Size (\d+)", segment)
            if not header_match:
                continue
                
            mat_size = int(header_match.group(1))
            blk_size = int(header_match.group(2))
            
            # 2. Find ALL CPU times (Tests 1-5)
            cpu_times = [float(m.group(1)) for m in rx_cpu.finditer(segment)]
            
            # 3. Find ALL GPU times
            gpu_matches = []
            for m in rx_gpu.finditer(segment):
                gpu_matches.append({
                    'type': m.group('type'),
                    'compute': float(m.group('compute')),
                    'total': float(m.group('total'))
                })
            
            # 4. Align and Flatten Data
            types = ['Unified', 'Shared', 'Global']
            
            for mem_type in types:
                relevant_matches = [m for m in gpu_matches if m['type'] == mem_type]
                
                for i, run in enumerate(relevant_matches):
                    row = {
                        'Matrix_Size': mat_size,
                        'Block_Size': blk_size,
                        'Memory_Type': mem_type,
                        'GPU_Compute': run['compute'],
                        'GPU_Total': run['total'],
                        'CPU_Time': cpu_times[i] if i < len(cpu_times) else None
                    }
                    data.append(row)
            
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return pd.DataFrame()

    return pd.DataFrame(data)

# ==========================================
# PART 2: PLOTTING & ANALYSIS
# ==========================================
def generate_analysis(df):
    if df.empty:
        print("No data parsed.")
        return

    # ---------------------------------------------------------
    # STEP 1: SCATTER PLOT & STATISTICS (Global vs Shared)
    # Configuration: Block 32, Matrix 4000
    # ---------------------------------------------------------
    
    # Create a copy for milliseconds conversion
    step1_df = df[(df['Block_Size'] == 32) & (df['Matrix_Size'] == 4000)].copy()
    if not step1_df.empty:
        step1_df['Compute_Time_ms'] = step1_df['GPU_Compute'] * 1000

        # Separate Data
        global_times = step1_df[step1_df["Memory_Type"] == "Global"]["Compute_Time_ms"]
        shared_times = step1_df[step1_df["Memory_Type"] == "Shared"]["Compute_Time_ms"]
        
        # 1. Print Statistics
        print("\n" + "="*50)
        print("STEP 1: STATISTICS (Block 32, Matrix 4000)")
        print("="*50)
        
        # Calculate Mean and Std Dev
        g_mean, g_std = global_times.mean(), global_times.std()
        s_mean, s_std = shared_times.mean(), shared_times.std()
        
        print(f"Global Memory:  {g_mean:.3f} ms ± {g_std:.3f} ms")
        print(f"Shared Memory:  {s_mean:.3f} ms ± {s_std:.3f} ms")
        print("="*50)

        # 2. Generate Plot
        plt.figure(figsize=(10, 2.5))

        # Plot Means (Large circles)
        plt.scatter(g_mean, 0, s=200, c='tab:blue', edgecolors="black", linewidths=2, zorder=3, label="Mean Global")
        plt.scatter(s_mean, 0, s=200, c='tab:orange', edgecolors="black", linewidths=2, zorder=3, label="Mean Shared")

        # Plot Individual Points (Smaller dots)
        jitter_g = np.random.normal(0, 0.02, size=len(global_times))
        jitter_s = np.random.normal(0, 0.02, size=len(shared_times))

        plt.scatter(global_times, jitter_g, alpha=0.6, s=100, c='tab:blue', label="Global Samples")
        plt.scatter(shared_times, jitter_s, alpha=0.6, s=100, c='tab:orange', label="Shared Samples")

        # Formatting
        plt.yticks([])  # Hide Y axis
        plt.xlabel(r"Compute Time $t$ (ms)", fontsize=12)
        plt.title("Step 1: Memory Performance Distribution (N=4000, B=32)")
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        
        # Adjust X-axis range
        all_vals = pd.concat([global_times, shared_times])
        margin = (all_vals.max() - all_vals.min()) * 0.5
        plt.xlim(all_vals.min() - margin, all_vals.max() + margin)
        
        plt.tight_layout()
        plt.savefig('measurements.png', dpi=300, bbox_inches='tight')
        print("Saved Plot: measurements.png")

    # ---------------------------------------------------------
    # STEP 2a: EXECUTION TIME TABLE (Global Memory)
    # ---------------------------------------------------------
    global_df = df[df['Memory_Type'] == 'Global']
    
    step2a_table = global_df.pivot_table(
        index='Matrix_Size', 
        columns='Block_Size', 
        values='GPU_Compute', 
        aggfunc='mean'
    )
    step2a_table.columns = [f"T_block_{c}" for c in step2a_table.columns]
    step2a_table.index.name = "N"
    
    print("\n" + "="*60)
    print("STEP 2a: Average Execution Time (Global Kernel) [seconds]")
    print("="*60)
    print(step2a_table)
    print(step2a_table.to_latex(float_format="%.6f"))

    # ---------------------------------------------------------
    # STEP 2b: EXECUTION TIME TABLE (Shared Memory)
    # ---------------------------------------------------------
    shared_df = df[df['Memory_Type'] == 'Shared']
    
    step2b_table = shared_df.pivot_table(
        index='Matrix_Size', 
        columns='Block_Size', 
        values='GPU_Compute', 
        aggfunc='mean'
    )
    step2b_table.columns = [f"T_block_{c}" for c in step2b_table.columns]
    step2b_table.index.name = "N"
    
    print("\n" + "="*60)
    print("STEP 2b: Average Execution Time (Shared Kernel) [seconds]")
    print("="*60)
    print(step2b_table)
    print(step2b_table.to_latex(float_format="%.6f"))

    # ---------------------------------------------------------
    # STEP 3: SPEEDUP ANALYSIS
    # ---------------------------------------------------------
    # Note: Using N=4000 as representative, though user mentioned 2000 in snippet. 
    # Sticking to N=4000 for high-load comparison.
    subset = df[(df['Matrix_Size'] == 4000) & (df['Block_Size'] == 32)]
    if not subset.empty:
        means = subset.groupby('Memory_Type')[['CPU_Time', 'GPU_Compute', 'GPU_Total']].mean()
        
        speedup_df = pd.DataFrame(index=means.index)
        speedup_df['Speedup_Compute'] = means['CPU_Time'] / means['GPU_Compute']
        speedup_df['Speedup_Total'] = means['CPU_Time'] / means['GPU_Total']
        
        print("\n" + "="*60)
        print("STEP 3: Speedup vs CPU (N=4000, Block=32)")
        print("="*60)
        print(speedup_df)

    # ---------------------------------------------------------
    # EXTRA PLOTS: Scaling & Speedup Evolution
    # ---------------------------------------------------------
    print("\n" + "="*60)
    print("GENERATING SCALING & SPEEDUP PLOTS...")
    print("="*60)

    # We filter for Block Size 32 as the baseline for scaling comparison
    scaling_subset = df[df['Block_Size'] == 32]
    
    # We aggregate means for each Matrix Size
    # This automatically handles whatever Matrix Sizes are present in the log (50...8000)
    agg_df = scaling_subset.groupby(['Matrix_Size', 'Memory_Type'])[['CPU_Time', 'GPU_Compute', 'GPU_Total']].mean().reset_index()

    # Pivot to get columns for each Memory Type
    pivot_compute = agg_df.pivot(index='Matrix_Size', columns='Memory_Type', values='GPU_Compute')
    pivot_cpu = agg_df.pivot(index='Matrix_Size', columns='Memory_Type', values='CPU_Time')
    # CPU time should be the same regardless of GPU memory type, so we just take one column (e.g., Global)
    # If CPU time varies slightly per run, averaging them is fine.
    cpu_times = pivot_cpu['Global'] 

    # Prepare data for plotting
    matrix_sizes = pivot_compute.index.values
    t_global = pivot_compute['Global'].values
    t_shared = pivot_compute['Shared'].values
    t_unified = pivot_compute['Unified'].values
    t_cpu = cpu_times.values

    # Plot Generation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- PLOT 1: Execution Time Scaling (Log-Log) ---
    ax1.plot(matrix_sizes, t_global, 'o-', label='Global Memory', linewidth=2)
    ax1.plot(matrix_sizes, t_shared, 's-', label='Shared Memory', linewidth=2)
    ax1.plot(matrix_sizes, t_unified, '^-', label='Unified Memory', linewidth=2) # Added Unified scaling
    ax1.plot(matrix_sizes, t_cpu, 'k--', label='CPU (Sequential)', linewidth=2, alpha=0.6)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Matrix Size ($N$)', fontsize=12)
    ax1.set_ylabel('Execution Time (s)', fontsize=12)
    ax1.set_title('Execution Time Scaling (Log-Log)', fontsize=14)
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend()

    # --- PLOT 2: Speedup vs Matrix Size ---
    # Calculate Speedups (Compute Only for this view to show kernel efficiency)
    speedup_global = t_cpu / t_global
    speedup_shared = t_cpu / t_shared
    speedup_unified = t_cpu / t_unified

    ax2.plot(matrix_sizes, speedup_global, 'o-', label='Global Speedup', linewidth=2)
    ax2.plot(matrix_sizes, speedup_shared, 's-', label='Shared Speedup', linewidth=2)
    ax2.plot(matrix_sizes, speedup_unified, '^-', label='Unified Speedup', linewidth=2)

    # Reference Line (Speedup = 1.0)
    ax2.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='CPU Parity')

    ax2.set_xscale('log')
    ax2.set_xlabel('Matrix Size ($N$)', fontsize=12)
    ax2.set_ylabel('Speedup ($T_{CPU} / T_{GPU}$)', fontsize=12)
    ax2.set_title('Compute Speedup Evolution', fontsize=14)
    ax2.grid(True, which="both", ls="-", alpha=0.2)
    ax2.legend()

    plt.tight_layout()
    plt.savefig('performance_analysis.png', dpi=300)
    print("Saved plots to performance_analysis.png")


# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    # Ensure matplotlib doesn't try to open a window if on a server
    plt.switch_backend('Agg') 
    
    df_robust = parse_lab_results('lab_results.out')
    generate_analysis(df_robust)
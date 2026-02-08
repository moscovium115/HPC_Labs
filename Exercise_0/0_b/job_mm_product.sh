#!/bin/bash
#SBATCH --job-name="MM_Scaling"
#SBATCH --time=00:15:00
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=1G
#SBATCH --account=Education-EEMCS-Courses-WI4049TU

# Get the names of the two nodes allocated to this job (e.g., node01, node02)
NODELIST=($(scontrol show hostnames $SLURM_JOB_NODELIST))
NODE_A=${NODELIST[0]}
NODE_B=${NODELIST[1]}

CUR_DIR_NAME=$(basename "$PWD")
DATA_DIR="${CUR_DIR_NAME}_data"

if [ ! -d "$DATA_DIR" ]; then
    mkdir "$DATA_DIR"
fi

echo "----------------------------------"
echo "Allocated Nodes: $NODE_A and $NODE_B"
echo "----------------------------------"

mpicc MM-parallel.c -o mm_parallel

echo "Starting Matrix Scaling Tests"
echo "--------------------------------------------"

P_VALUES="1 2 8 24 40 41 48 64"

for P in $P_VALUES
do
    echo "Processing P=$P..."

    # --- CASE 1: Shared Memory (Fit on 1 Node) ---
    if [ "$P" -le 40 ]; then
        # Just use standard block distribution on 1 node.
        # It will put all P tasks on NODE_A.
        echo "  -> Mode: Shared Memory (Single Node)"
        
        srun --nodes=1 --ntasks=$P --distribution=block ./mm_parallel > "$DATA_DIR/mm_out_P${P}.txt"

    # --- CASE 2: Distributed Memory (Force Saturation) ---
    else
        echo "  -> Mode: Distributed (Saturating Node A)"
        
        # Create a custom hostfile for 'arbitrary' distribution
        # We want: 40 tasks on NODE_A, remainder on NODE_B
        HOSTFILE="hostfile_P${P}.txt"
        > "$HOSTFILE" # Create/Clear file
        
        # Add Node A 40 times
        for ((i=1; i<=40; i++)); do
            echo "$NODE_A" >> "$HOSTFILE"
        done
        
        # Add Node B (P-40) times
        REMAINDER=$((P - 40))
        for ((i=1; i<=REMAINDER; i++)); do
            echo "$NODE_B" >> "$HOSTFILE"
        done
        
        # Tell SLURM to use this specific file
        export SLURM_HOSTFILE=./$HOSTFILE
        
        # Layout Verification
        # We use --distribution=arbitrary so it follows our file exactly
        srun --nodes=2 --ntasks=$P --distribution=arbitrary hostname | sort | uniq -c
        
        # Run the actual code
        srun --nodes=2 --ntasks=$P --distribution=arbitrary ./mm_parallel > "$DATA_DIR/mm_out_P${P}.txt"
        
        # Cleanup temp file
        rm "$HOSTFILE"
    fi
done

# --- Combining Results ---

FINAL_OUTPUT="${DATA_DIR}/mm_out_all.txt"
echo "Consolidating all results into $FINAL_OUTPUT..."

> "$FINAL_OUTPUT"

for P in $P_VALUES
do
    echo "========================================" >> "$FINAL_OUTPUT"
    echo " RESULTS FOR P = $P "                     >> "$FINAL_OUTPUT"
    echo "========================================" >> "$FINAL_OUTPUT"
    
    if [ -f "$DATA_DIR/mm_out_P${P}.txt" ]; then
        cat "$DATA_DIR/mm_out_P${P}.txt" >> "$FINAL_OUTPUT"
    else
        echo "Error: Output file for P=$P not found." >> "$FINAL_OUTPUT"
    fi
    
    echo -e "\n" >> "$FINAL_OUTPUT"
done

echo "Tests complete. All data saved to $FINAL_OUTPUT"
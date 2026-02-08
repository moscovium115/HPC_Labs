#!/bin/sh
#SBATCH --job-name=Poisson_Ex1.2.2
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# --- 1. CONFIGURATION ---
# Create the 'data' directory if it doesn't exist
# This is safe; it won't error if the directory already exists
mkdir -p data

# Configure input.dat for Exercise 1.2.2 (Gridsize 100, 4 processors) 
# Modifying in place, assuming file exists.
echo "Configuring input.dat (100x100, max_iter=20000)..."
sed -i 's/^nx:.*/nx: 100/' input.dat
sed -i 's/^ny:.*/ny: 100/' input.dat
# Increase max iterations to prevent premature stopping
sed -i 's/^max iterations:.*/max iterations: 20000/' input.dat

echo "======================================================================="
echo "Exercise 1.2.2: Optimization of Omega"
echo "Topology: 4x1 | Grid: 100x100 | Processes: 4"
echo "Log files saved to: ./data/"
echo "======================================================================="
printf "%-10s | %-15s | %-15s | %-15s\n" "Omega" "Iterations" "Time (s)" "Efficiency (%)"
echo "-----------------------------------------------------------------------"

# --- 2. LOOP THROUGH OMEGA VALUES ---
# [cite_start]Strategically chosen values between 1.90 and 1.99 (5 runs maximal) [cite: 36]
for W in 1.00 1.95
do
    # Modify the C source code directly to set the new omega
    sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = $W;/" MPI_Poisson_1_2_1.c

    # Recompile
    mpicc MPI_Poisson_1_2_1.c -o Poisson.exe -lm

    # Run on 4 processors with 4x1 topology 
    # Redirect output to the 'data' folder
    srun -n 4 ./Poisson.exe 4 1 > data/output_omega_${W}.log
    
    # --- 3. DATA EXTRACTION ---
    # Note: We must now grep from the file inside 'data/'
    
    # Get Iterations
    ITER=$(grep "Number of iterations" data/output_omega_${W}.log | head -n 1 | awk '{print $NF}')
    
    # Get Wallclock Time
    TIME=$(grep "Elapsed Wtime" data/output_omega_${W}.log | head -n 1 | awk '{print $4}')
    
    # Get CPU Efficiency
    EFF=$(grep "Elapsed Wtime" data/output_omega_${W}.log | head -n 1 | grep -o "[0-9.]*%" | tr -d '%')

    # Print formatted row
    printf "%-10s | %-15s | %-15s | %-15s\n" "$W" "$ITER" "$TIME" "$EFF"
done

# --- 4. RESET ---
echo "-----------------------------------------------------------------------"
echo "Resetting omega to 1.0 (Standard Gauss-Seidel)..."
# Reset the C file back to default omega = 1.0
sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = 1.0;/" MPI_Poisson_1_2_1.c

echo "Done."
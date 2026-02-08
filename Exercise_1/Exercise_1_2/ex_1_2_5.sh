#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.5
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --nodes=1              
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# --- SETUP ---
mkdir -p data_ex_1_2_5
C_FILE="MPI_Poisson_1_2_5.c"
EXE="Poisson_Convergence.exe"
DATA_FILE="data_ex_1_2_5/convergence_results.csv"


# Ensure Omega is set (using 1.95 as a good baseline)
sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = 1.95;/" $C_FILE

# Compile
echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- 2. EXPERIMENT LOOP ---
echo "Grid,Iterations,Time" > $DATA_FILE

echo "======================================================================="
echo "Running Convergence Experiments (Fixed Omega=1.95)"
echo "Output CSV: $DATA_FILE"
echo "======================================================================="

# Using 4x1 Topology as our standard reference (since it was the fastest)
PX=4
PY=1

# Loop Grid Sizes: 200 up to 2000
for GRID in 200 400 800 1000 1500
do
    echo "--> Running Grid: ${GRID}x${GRID}..."
    
    # --- MODIFY INPUT.DAT ---
    # 1. Update Grid Size
    sed -i "s/^nx:.*/nx: $GRID/" input.dat
    sed -i "s/^ny:.*/ny: $GRID/" input.dat
    
    # 2. RESTORE Precision Goal (Crucial for 1.2.5)
    sed -i "s/^precision goal:.*/precision goal: 1e-4/" input.dat
    
    # 3. Set Max Iterations high enough to allow convergence
    sed -i "s/^max iterations:.*/max iterations: 50000/" input.dat

    # --- RUN EXPERIMENT ---
    LOGFILE="data_ex_1_2_5/grid${GRID}.log"
    srun -n 4 ./$EXE $PX $PY > $LOGFILE

    # --- PARSE OUTPUT ---
    # Extract iteration count. 
    # Output format: "(0) Number of iterations : 1234"
    # We take the max iteration count reported (should be same for all ranks)
    ITERS=$(grep "Number of iterations" $LOGFILE | awk '{print $6}' | sort -nr | head -n 1)
    
    # Extract Total Time (optional, but interesting)
    # Output format: "(0) Elapsed Wtime ..."
    TIME=$(grep "Elapsed Wtime" $LOGFILE | awk '{print $4}' | sort -nr | head -n 1)

    echo "${GRID},${ITERS},${TIME}" >> $DATA_FILE
    echo "    Converged in $ITERS iterations ($TIME s)"
done

echo "======================================================================="
echo "Done."

# --- CLEANUP ---
# Restore C Code to original state (optional)
# sed -i "s/\/\/global_delta=1;/global_delta=1;/" $C_FILE
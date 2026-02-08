#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.9
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
mkdir -p data_ex_1_2_9
C_FILE="MPI_Poisson_1_2_9.c"
EXE="Poisson_Optimized.exe"
DATA_FILE="data_ex_1_2_9/optimization_result.csv"

# --- COMPILE ---
if [ ! -f "$C_FILE" ]; then
    echo "Error: $C_FILE not found!"
    exit 1
fi

echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- RUN ---
echo "Version,Iterations,Time" > $DATA_FILE

echo "Running Optimized Code (800x800)..."
LOGFILE="data_ex_1_2_9/run.log"

# Run with 4 processors (4x1 Topology)
# Assumes input.dat is present in the current directory
srun -n 4 ./$EXE 4 1 > $LOGFILE

# Parse Results
ITERS=$(grep "Number of iterations" $LOGFILE | awk '{print $6}' | sort -nr | head -n 1)
TIME=$(grep "Elapsed Wtime" $LOGFILE | awk '{print $4}' | sort -nr | head -n 1)

echo "Optimized,$ITERS,$TIME" >> $DATA_FILE

echo "Done. Result: $ITERS iters in $TIME s"
echo "Results saved to $DATA_FILE"
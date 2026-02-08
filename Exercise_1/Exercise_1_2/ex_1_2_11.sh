#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.11
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:10:00
#SBATCH --nodes=1              
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# --- SETUP ---
mkdir -p data_ex_1_2_11
C_FILE="MPI_Poisson_1_2_11.c"
EXE="Poisson_Comm.exe"
DATA_FILE="data_ex_1_2_11/comm_analysis.csv"

# --- COMPILE ---
if [ ! -f "$C_FILE" ]; then
    echo "Error: $C_FILE not found!"
    exit 1
fi

echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- RUN EXPERIMENTS ---
echo "Topology,GridSize,AvgCommTime" > $DATA_FILE
SIZES="200 400 800 1000"

echo "Running Communication Analysis..."

# --- 1. Run 4x1 Topology ---
for N in $SIZES; do
    echo "  -> 4x1 with Grid $N..."
    echo "nx: $N" > input.dat
    echo "ny: $N" >> input.dat
    echo "precision goal: 1e-4" >> input.dat 
    echo "max iterations: 1000" >> input.dat 
    echo "source: 0.35 0.70 4.0" >> input.dat
    echo "source: 0.625 0.75 4.0" >> input.dat
    echo "source: 0.375 0.25 4.0" >> input.dat

    # Run and capture output
    OUTPUT=$(srun -n 4 ./$EXE 4 1)
    
    # FIX: Use $12 to grab the number AFTER "COMM_TIME"
    TIME=$(echo "$OUTPUT" | grep "RANK 0" | awk '{print $12}')
    
    echo "4x1,$N,$TIME" >> $DATA_FILE
done

# --- 2. Run 2x2 Topology ---
for N in $SIZES; do
    echo "  -> 2x2 with Grid $N..."
    echo "nx: $N" > input.dat
    echo "ny: $N" >> input.dat
    echo "precision goal: 1e-4" >> input.dat 
    echo "max iterations: 1000" >> input.dat 
    echo "source: 0.35 0.70 4.0" >> input.dat 
    echo "source: 0.625 0.75 4.0" >> input.dat
    echo "source: 0.375 0.25 4.0" >> input.dat

    OUTPUT=$(srun -n 4 ./$EXE 2 2)
    # FIX: Use $12 to grab the number AFTER "COMM_TIME"
    TIME=$(echo "$OUTPUT" | grep "RANK 0" | awk '{print $12}')
    
    echo "2x2,$N,$TIME" >> $DATA_FILE
done

# --- RESTORE ORIGINAL INPUT.DAT ---
echo "Restoring original input.dat..."
echo "nx: 100" > input.dat
echo "ny: 100" >> input.dat
echo "precision goal: 0.0001" >> input.dat
echo "max iterations: 5000" >> input.dat
echo "source: 0.35 0.70 4.0" >> input.dat
echo "source: 0.625 0.75 4.0" >> input.dat
echo "source: 0.375 0.25 4.0" >> input.dat

echo "Done. Results saved to $DATA_FILE"
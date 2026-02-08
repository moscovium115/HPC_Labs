#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.7
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
mkdir -p data_ex_1_2_7
C_FILE="MPI_Poisson_1_2_7.c"
EXE="Poisson_ErrorFreq.exe"
DATA_FILE="data_ex_1_2_7/freq_results.csv"

# --- PREPARE INPUT.DAT ---
echo "Configuring input.dat for benchmark (800x800)..."
sed -i 's/^nx:.*/nx: 800/' input.dat
sed -i 's/^ny:.*/ny: 800/' input.dat
sed -i 's/^max iterations:.*/max iterations: 50000/' input.dat
sed -i 's/^precision goal:.*/precision goal: 1e-4/' input.dat

echo "Frequency,Iterations,Time" > $DATA_FILE

# --- 1. RUN STANDARD FREQUENCIES (1, 10, 100) ---
for FREQ in 1 10 100
do
    echo "------------------------------------------"
    echo "Benchmarking Error Check Frequency: $FREQ"

    # Modify C code
    sed -i "s/int error_check_freq=[0-9]*;/int error_check_freq=$FREQ;/" $C_FILE
    mpicc $C_FILE -o $EXE -lm

    # Run
    OUTPUT=$(srun -n 4 ./$EXE 4 1)
    
    # Parse
    ITERS=$(echo "$OUTPUT" | grep "Number of iterations" | head -n 1 | awk '{print $6}')
    TIME=$(echo "$OUTPUT" | grep "Elapsed Wtime" | sort -nr | head -n 1 | awk '{print $4}')
    
    echo "$FREQ,$ITERS,$TIME" >> $DATA_FILE
    echo "  -> Result: $ITERS iters in $TIME s"
done

# --- 2. RUN "ZERO CHECK" CASE ---
echo "------------------------------------------"
echo "Benchmarking ZERO Error Checks (Fixed 3003 Iters)"

# Set Frequency to 10000 (so it never checks)
sed -i "s/int error_check_freq=[0-9]*;/int error_check_freq=10000;/" $C_FILE
mpicc $C_FILE -o $EXE -lm

# Force input.dat to stop exactly at 3003
sed -i 's/^max iterations:.*/max iterations: 3003/' input.dat

# Run
OUTPUT=$(srun -n 4 ./$EXE 4 1)

ITERS=$(echo "$OUTPUT" | grep "Number of iterations" | head -n 1 | awk '{print $6}')
TIME=$(echo "$OUTPUT" | grep "Elapsed Wtime" | sort -nr | head -n 1 | awk '{print $4}')

# Label this as "0" or "Inf" in the CSV
echo "0,$ITERS,$TIME" >> $DATA_FILE
echo "  -> Result: $ITERS iters in $TIME s"

# --- CLEANUP ---
# Restore defaults
sed -i "s/int error_check_freq=[0-9]*;/int error_check_freq=1;/" $C_FILE
sed -i 's/^nx:.*/nx: 100/' input.dat
sed -i 's/^ny:.*/ny: 100/' input.dat
sed -i 's/^max iterations:.*/max iterations: 5000/' input.dat

echo "------------------------------------------"
echo "Done. Results saved to $DATA_FILE"
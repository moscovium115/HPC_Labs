#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.8
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:20:00
#SBATCH --nodes=1              
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# --- SETUP ---
mkdir -p data_ex_1_2_8
C_FILE="MPI_Poisson_1_2_8.c"
EXE="Poisson_Sweeps.exe"
DATA_FILE="data_ex_1_2_8/sweep_results.csv"

# --- 1. COMPILE ---
if [ ! -f "$C_FILE" ]; then
    echo "Error: $C_FILE not found!"
    exit 1
fi

echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- 2. RUN EXPERIMENTS ---
echo "Sweeps,Iterations,Time" > $DATA_FILE

echo "======================================================================="
echo "Running Sweep Optimization (k=1..10) on 800x800 Grid"
echo "Output CSV: $DATA_FILE"
echo "======================================================================="

# Loop through Sweeps k = 1 to 10
# We use the 4x1 topology (Arguments: 4 1) which was the fastest
for K in {1..10}
do
    echo "--> Running with Sweeps per Exchange = $K..."
    
    LOGFILE="data_ex_1_2_8/sweeps_$K.log"
    
    # Execution arguments: [Px] [Py] [Sweeps]
    srun -n 4 ./$EXE 4 1 $K > $LOGFILE
    
    # Parse Results
    # Extract iteration count (The line format is "(rank) Number of iterations : ...")
    ITERS=$(grep "Number of iterations" $LOGFILE | awk '{print $6}' | sort -nr | head -n 1)
    
    # Extract Time
    TIME=$(grep "Elapsed Wtime" $LOGFILE | awk '{print $4}' | sort -nr | head -n 1)
    
    echo "$K,$ITERS,$TIME" >> $DATA_FILE
    echo "    Result: $ITERS iterations in $TIME s"
done

echo "======================================================================="
echo "Done."
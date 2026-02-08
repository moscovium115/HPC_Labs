#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.3_Trace
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
# Create directory for logs and results
mkdir -p data_ex_1_2_3

C_FILE="MPI_Poisson_1_2_3.c"
EXE="Poisson_1_2_3.exe"
DATA_FILE="data_ex_1_2_3/scaling_trace.csv"

# 1. Hardcode Omega = 1.95 in the C file
sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = 1.95;/" $C_FILE

# 2. Compile
echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# 3. Initialize CSV File
echo "Topology,Grid,Iterations,Time" > $DATA_FILE

echo "======================================================================="
echo "Running Trace Experiments (Single Run per Config)"
echo "Output CSV: $DATA_FILE"
echo "======================================================================="

# --- EXPERIMENT LOOPS ---

# Loop 1: Topologies (Px Py)
for TOPOLOGY in "4 1" "2 2" "1 4"
do
    set -- $TOPOLOGY
    PX=$1
    PY=$2

    # Loop 2: Grid Sizes
    # Using larger max iterations (2000) to get 20 data points (printing every 100)
    for GRID in 200 400 800 1200 1500
    do
        echo "--> Config: ${PX}x${PY} | Grid: ${GRID}x${GRID}"
        
        # --- MODIFY INPUT.DAT SAFELY ---
        # 1. Update Grid Size
        sed -i "s/^nx:.*/nx: $GRID/" input.dat
        sed -i "s/^ny:.*/ny: $GRID/" input.dat
        
        # 2. Set Precision to 0.0 (Ensure it runs the full duration)
        sed -i "s/^precision goal:.*/precision goal: 0.0/" input.dat
        
        # 3. Set Max Iterations to 2000 (To generate enough 'MEASURE' points)
        sed -i "s/^max iterations:.*/max iterations: 2000/" input.dat

        # --- RUN EXPERIMENT ---
        LOGFILE="data_ex_1_2_3/topo${PX}x${PY}_grid${GRID}.log"
        
        srun -n 4 ./$EXE $PX $PY > $LOGFILE

        # --- PARSE OUTPUT ---
        # Your code prints: "MEASURE: Iterations <count> Time <time>"
        # We grep this line and extract the numbers to the CSV
        grep "MEASURE:" $LOGFILE | while read -r LABEL LABEL2 ITER LABEL3 TIME
        do
            # Format: Topology,Grid,Iterations,Time
            echo "${PX}x${PY},${GRID},${ITER},${TIME}" >> $DATA_FILE
        done
    done
done

echo "======================================================================="
echo "Experiments Complete. Restoring defaults."

# --- CLEANUP ---

# Restore Omega to 1.0
sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = 1.0;/" $C_FILE

# Restore input.dat defaults
sed -i "s/^nx:.*/nx: 100/" input.dat
sed -i "s/^ny:.*/ny: 100/" input.dat
sed -i "s/^precision goal:.*/precision goal: 1e-4/" input.dat
sed -i "s/^max iterations:.*/max iterations: 20000/" input.dat

# Recompile clean
mpicc $C_FILE -o $EXE -lm
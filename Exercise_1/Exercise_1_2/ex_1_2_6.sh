#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.6
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
mkdir -p data_ex_1_2_6
C_FILE="MPI_Poisson_1_2_6.c"
EXE="Poisson_Error.exe"
DATA_FILE="data_ex_1_2_6/error_trace_800.csv"

# --- 1. PATCH C CODE (Add Error Logging) ---
# We inject the printf statement into the loop if it's not there
# This sed looks for the MPI_Allreduce line and adds the print after it
# (Ideally, just use the manually modified C file above, but this is a quick patch)
# BETTER OPTION: Just manually edit the C file as shown above before running this script.

# Compile
echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- 2. CONFIGURE INPUT ---
# 800x800 Grid, Omega=1.95, Tolerance=1e-4
echo "nx: 800" > input.dat
echo "ny: 800" >> input.dat
echo "precision goal: 1e-4" >> input.dat
echo "max iterations: 50000" >> input.dat
echo "source: 0.35 0.70 4.0" >> input.dat
echo "source: 0.625 0.75 4.0" >> input.dat
echo "source: 0.375 0.25 4.0" >> input.dat

# Ensure Omega is 1.95 in the C file (if hardcoded) or add it to input if supported
sed -i "s/double[[:space:]]*omega[[:space:]]*=[[:space:]]*[0-9.]*;/double omega = 1.95;/" $C_FILE
mpicc $C_FILE -o $EXE -lm

# --- 3. RUN EXPERIMENT ---
echo "Iteration,Error" > $DATA_FILE
echo "Running 800x800 Error Monitor..."

srun -n 4 ./$EXE 4 1 > data_ex_1_2_6/run.log

# Extract the ERROR_LOG lines
grep "ERROR_LOG" data_ex_1_2_6/run.log | awk '{print $2 "," $3}' >> $DATA_FILE

echo "Done. Data saved to $DATA_FILE."
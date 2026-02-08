#!/bin/sh
#SBATCH --job-name=Poisson_Step9
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:01:00
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# Compile the code
mpicc MPI_Poisson.c -o Poisson.exe -lm

echo "----------------------------------------"
echo "BASELINE: 1 Process (1x1 grid)"
# Run on 1 processor to get the sequential iteration count
srun -n 1 ./Poisson.exe 1 1 > output_baseline.log
grep "Number of iterations" output_baseline.log

echo "----------------------------------------"
echo "TEST: 4 Processes (2x2 grid)"
# Run on 4 processors. If Step 9 is correct, this will NOT hang 
# and will match the baseline iterations exactly.
srun -n 4 ./Poisson.exe 2 2 > output_parallel.log
grep "Number of iterations" output_parallel.log
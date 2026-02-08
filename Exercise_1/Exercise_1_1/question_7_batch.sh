#!/bin/sh
#SBATCH --job-name=Poisson_Step7
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --ntasks=4             # Max tasks we will need (for 2x2)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# Compile
mpicc MPI_Poisson.c -o Poisson.exe -lm

echo "----------------------------------------"
echo "RUN 1: 2 Processes (2x1 grid)"
# Split domain into 2 vertical strips
srun -n 2 ./Poisson.exe 2 1 > output_2.log
grep "Number of iterations" output_2.log

echo "----------------------------------------"
echo "RUN 2: 3 Processes (3x1 grid)"
# Split domain into 3 vertical strips (Left, Middle, Right)
srun -n 3 ./Poisson.exe 3 1 > output_3.log
grep "Number of iterations" output_3.log

echo "----------------------------------------"
echo "RUN 3: 4 Processes (2x2 grid)"
# Split domain into 4 quadrants
srun -n 4 ./Poisson.exe 2 2 > output_4.log
grep "Number of iterations" output_4.log
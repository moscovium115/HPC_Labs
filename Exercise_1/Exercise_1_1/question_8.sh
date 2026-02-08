#!/bin/sh
#SBATCH --job-name=Poisson_Step7
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:01:00
#SBATCH --ntasks=4             # Max tasks we will need (for 2x2)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# Compile
mpicc MPI_Poisson.c -o Poisson.exe -lm



echo "----------------------------------------"
echo "RUN 3: 4 Processes (2x2 grid)"
# Split domain into 4 quadrants
srun -n 4 ./Poisson.exe 2 2 > output_4.log
grep "Number of iterations" output_4.log
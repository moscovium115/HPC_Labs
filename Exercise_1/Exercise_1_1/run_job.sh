#!/bin/sh
#SBATCH --job-name=Poisson_Step1
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# Compile the program
mpicc MPI_Poisson.c -o Poisson.exe -lm

# Run on 2 processors
# We output to a log file to verify the prints from both clones
srun ./Poisson.exe  1 > output.log
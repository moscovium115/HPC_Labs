#!/bin/bash
#SBATCH --job-name=FEM_Ex4_5
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

# Define output filename
OUTPUT_FILE="ex4_6_output.txt"

# Clear the old output file
> "$OUTPUT_FILE"

# Recompile to be safe
echo "Compiling..."
make

# Loop through the three sizes with the 'adapt' keyword
for SIZE in 100 200 400; do
    echo "----------------------------------------" | tee -a "$OUTPUT_FILE"
    echo "Running ADAPTIVE simulation for Grid: $SIZE x $SIZE" | tee -a "$OUTPUT_FILE"
    echo "----------------------------------------" | tee -a "$OUTPUT_FILE"

    # 1. Generate the ADAPTIVE grid
    # Usage: ./GridDist <Px> <Py> <Nx> <Ny> adapt
    echo "Generating adaptive grid for $SIZE x $SIZE..."
    ./GridDist 2 2 $SIZE $SIZE adapt

    # 2. Run the solver
    echo "Running solver..."
    mpirun -n 4 ./MPI_Fempois >> "$OUTPUT_FILE"
    
    echo "" >> "$OUTPUT_FILE"
done

echo "Done! Adaptive results saved to $OUTPUT_FILE"
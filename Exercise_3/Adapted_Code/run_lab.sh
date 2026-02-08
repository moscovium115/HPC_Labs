#!/bin/sh
#SBATCH --job-name="PowerGPU_Lab3"
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --partition=gpu-a100-small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=lab_results.out

# 1. Clean and Load Environment
module purge
module load 2025
module load nvhpc

# 2. Fix for "cannot find crtbegin.o"
export LIBRARY_PATH=/usr/lib/gcc/x86_64-redhat-linux/8/:$LIBRARY_PATH

echo "=========================================="
echo "STARTING EXPERIMENTS WITH RECOMPILATION"
echo "=========================================="

# Loop through Block Sizes first (so we recompile fewer times)
for BS in 32 64 100
do
    echo ""
    echo "##########################################"
    echo "SETTING CONSTANT BLOCK_SIZE = $BS"
    echo "##########################################"

    # --- THE TRICK: MODIFY C++ CODE AND RECOMPILE ---
    
    # 1. Use sed to replace "const int BLOCK_SIZE = ...;" with new value
    # We look for the line starting with "const int BLOCK_SIZE =" and replace the number
    sed -i "s/const int BLOCK_SIZE =.*;/const int BLOCK_SIZE =$BS;/g" power_gpu.cu

    # 2. Recompile
    echo "Recompiling for Block Size $BS..."
    nvcc -ccbin /usr/bin/g++ -O3 power_gpu.cu -o power_gpu

    # Safety Check
    if [ ! -f ./power_gpu ]; then
        echo "ERROR: Compilation failed for Block Size $BS. Exiting."
        # Reset file before exiting
        sed -i "s/const int BLOCK_SIZE =.*;/const int BLOCK_SIZE =32;/g" power_gpu.cu
        exit 1
    fi

    # 3. Run Experiments for this Block Size
    for N in 50 500 2000 4000 8000
    do
        echo ""
        echo ">>> RUNNING: Matrix Size $N | Block Size $BS"
        echo "------------------------------------------"
        # We pass --block_size too so the Host code knows how many threads to launch
        ./power_gpu --size $N --block_size $BS
        echo "------------------------------------------"
    done
done

# --- CLEANUP ---
echo ""
echo "Resetting code back to default..."
sed -i "s/const int BLOCK_SIZE =.*;/const int BLOCK_SIZE =32;/g" power_gpu.cu

echo "=========================================="
echo "ALL EXPERIMENTS FINISHED"
echo "=========================================="
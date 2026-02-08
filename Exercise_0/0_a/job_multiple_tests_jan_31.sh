#!/bin/bash
#SBATCH --job-name="PingPong_Both"
#SBATCH --time=00:05:00
#SBATCH --partition=compute

# Allocate 2 nodes by requesting 4 tasks with max 2 per node
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=Education-EEMCS-Courses-WI4049TU

CUR_DIR_NAME=$(basename "$PWD")
DATA_DIR="${CUR_DIR_NAME}_data"

if [ ! -d "$DATA_DIR" ]; then
    echo "Creating data directory: $DATA_DIR"
    mkdir "$DATA_DIR"
else
    echo "Data directory already exists: $DATA_DIR"
fi

echo "----------------------------------"

# Compile once
mpicc pingPong.c -o pingPong

echo "Allocated nodes (Total pool):"
scontrol show hostnames $SLURM_NODELIST
echo "-----------------------------"

########################################
# Test 2: DIFFERENT NODES (1 rank per node)
########################################
echo ">>> STARTING TEST 2: Different Nodes (Inter-node)"
echo "Verifying node distribution (Task ID : Hostname):"

# Added --label to see Rank ID
# Added sort to order the output
# Added tee to print to screen AND save to file
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 --label hostname | sort | tee "$DATA_DIR/ping_diff_hosts.txt"

echo "Running PingPong executable..."
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 ./pingPong > "$DATA_DIR/ping_diff.txt"
echo "Done. Output saved to $DATA_DIR/ping_diff.txt"
echo "-----------------------------"

########################################
# Test 1: SAME NODE (2 ranks on 1 node)
########################################
echo ">>> STARTING TEST 1: Same Node (Intra-node)"
echo "Verifying node distribution (Task ID : Hostname):"

srun --nodes=1 --ntasks=2 --ntasks-per-node=2 --label hostname | sort | tee "$DATA_DIR/ping_same_hosts.txt"

echo "Running PingPong executable..."
srun --nodes=1 --ntasks=2 --ntasks-per-node=2 ./pingPong > "$DATA_DIR/ping_same.txt"
echo "Done. Output saved to $DATA_DIR/ping_same.txt"

echo "-----------------------------"




echo "All tests finished."
#!/bin/bash
#SBATCH --job-name="PingPong_Combined"
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

# 1. COMPILE BOTH PROGRAMS
echo "Compiling Standard PingPong..."
mpicc pingPong.c -o pingPong

echo "Compiling SendRecv PingPong..."
mpicc pingPong_SendRecv.c -o pingPongExtra

echo "----------------------------------"
echo "Allocated nodes (Total pool):"
scontrol show hostnames $SLURM_NODELIST
echo "----------------------------------"


########################################
# PART 1: STANDARD PING PONG TESTS
########################################

echo "=== PART 1: STANDARD PING PONG ==="

# --- Test 1.1: DIFFERENT NODES ---
echo ">>> [Standard] Testing Different Nodes (Inter-node)"
# Verify nodes on screen and save to file
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 --label hostname | sort | tee "$DATA_DIR/ping_std_diff_hosts.txt"

echo "Running Standard PingPong (Diff)..."
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 ./pingPong > "$DATA_DIR/ping_std_diff.txt"

# --- Test 1.2: SAME NODE ---
echo ">>> [Standard] Testing Same Node (Intra-node)"
# Verify nodes on screen and save to file
srun --nodes=1 --ntasks=2 --ntasks-per-node=2 --label hostname | sort | tee "$DATA_DIR/ping_std_same_hosts.txt"

echo "Running Standard PingPong (Same)..."
srun --nodes=1 --ntasks=2 --ntasks-per-node=2 ./pingPong > "$DATA_DIR/ping_std_same.txt"


echo "----------------------------------"


########################################
# PART 2: SEND/RECV EXTRA TESTS
########################################

echo "=== PART 2: SEND/RECV PING PONG ==="

# --- Test 2.1: DIFFERENT NODES ---
echo ">>> [SendRecv] Testing Different Nodes (Inter-node)"
# Verify nodes on screen and save to file
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 --label hostname | sort | tee "$DATA_DIR/ping_diff_hosts_extra.txt"

echo "Running SendRecv PingPong (Diff)..."
srun --nodes=2 --ntasks=2 --ntasks-per-node=1 ./pingPongExtra > "$DATA_DIR/ping_diff_extra.txt"

# --- Test 2.2: SAME NODE ---
echo ">>> [SendRecv] Testing Same Node (Intra-node)"
# Verify nodes on screen and save to file
srun --nodes=1 --ntasks=2 --ntasks-per-node=2 --label hostname | sort | tee "$DATA_DIR/ping_same_hosts_extra.txt"

echo "Running SendRecv PingPong (Same)..."
srun --nodes=1 --ntasks=2 --ntasks-per-node=2 ./pingPongExtra > "$DATA_DIR/ping_same_extra.txt"

echo "----------------------------------"
echo "All standard and extra tests finished."
echo "Output files located in: $DATA_DIR"
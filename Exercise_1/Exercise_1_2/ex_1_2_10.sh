#!/bin/bash
#SBATCH --job-name=Poisson_Ex1.2.10
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:20:00
#SBATCH --nodes=1              
#SBATCH --ntasks=36            
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

# --- SETUP ---
mkdir -p data_ex_1_2_10
C_FILE="MPI_Poisson_1_2_10.c"
EXE="Poisson_Ratio.exe"
DATA_FILE="data_ex_1_2_10/ratio_results.csv"

# --- COMPILE ---
if [ ! -f "$C_FILE" ]; then
    echo "Error: $C_FILE not found!"
    exit 1
fi
echo "Compiling $C_FILE..."
mpicc $C_FILE -o $EXE -lm

# --- RUN EXPERIMENTS ---
echo "GridSize,Processes,CalcTime,CommTime,Ratio" > $DATA_FILE

GRID_SIZES="50 100 200 400 800"
PROCESS_COUNTS="1 4 9 16 25 36"

echo "Running Ratio Analysis (Strong Scaling Sweep)..."

for N in $GRID_SIZES
do
    echo "------------------------------------------"
    echo "Testing Grid Size: $N x $N"
    
    # Configure input.dat for this Grid Size
    echo "nx: $N" > input.dat
    echo "ny: $N" >> input.dat
    echo "precision goal: 1e-4" >> input.dat
    echo "max iterations: 1000" >> input.dat 
    echo "source: 0.35 0.70 4.0" >> input.dat
    echo "source: 0.625 0.75 4.0" >> input.dat
    echo "source: 0.375 0.25 4.0" >> input.dat

    for P in $PROCESS_COUNTS
    do
        # Define Square Topology (Px x Py)
        case $P in
            1)  PX=1; PY=1 ;;
            4)  PX=2; PY=2 ;;
            9)  PX=3; PY=3 ;;
            16) PX=4; PY=4 ;;
            25) PX=5; PY=5 ;;
            36) PX=6; PY=6 ;;
        esac

        # Run
        OUTPUT=$(srun -n $P ./$EXE $PX $PY)
        
        # Parse Output
        LINE=$(echo "$OUTPUT" | grep "^RESULT")
        CALC=$(echo "$LINE" | awk '{print $3}')
        COMM=$(echo "$LINE" | awk '{print $4}')
        
        # Calculate Ratio
        RATIO=$(python3 -c "print(f'{($COMM / $CALC):.4f}' if float('$CALC') > 0 else 0)")

        echo "$N,$P,$CALC,$COMM,$RATIO" >> $DATA_FILE
        echo "  -> P=$P: Calc=$CALC s | Comm=$COMM s | Ratio=$RATIO"
    done
done

echo "Done. Results saved to $DATA_FILE"
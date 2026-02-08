#!/bin/bash
#SBATCH --job-name=FEM_Ex4_4
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --nodes=2
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

set -euo pipefail

# --- output dirs ---
BASE_OUT="data_ex_4_4"
mkdir -p "$BASE_OUT"

NEIGH_OUT="${BASE_OUT}/neighbors_2x2_and_3x3.txt"
LOG_OUT="${BASE_OUT}/run_log.txt"
: > "$NEIGH_OUT"
: > "$LOG_OUT"

# Files GridDist expects in CWD
REQ_FILES=(sources.dat input.dat)

echo "Compiling..." | tee -a "$LOG_OUT"
make | tee -a "$LOG_OUT"
echo "" | tee -a "$LOG_OUT"

SIZE=200

stage_required_files () {
  local DIR=$1
  for f in "${REQ_FILES[@]}"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: required file '$f' not found in project root ($(pwd))" | tee -a "$LOG_OUT"
      exit 1
    fi
    cp -f "$f" "$DIR/$f"
  done
}

run_case () {
  local PX=$1
  local PY=$2
  local NP=$((PX*PY))
  local CASE_DIR="${BASE_OUT}/case_${PX}x${PY}_N${SIZE}"

  mkdir -p "$CASE_DIR"

  # Put sources.dat + input.dat into the case dir for GridDist
  stage_required_files "$CASE_DIR"

  echo "==================================================" >> "$NEIGH_OUT"
  echo "Neighbour counts for process grid ${PX}x${PY} (P=${NP}), grid ${SIZE}x${SIZE}" >> "$NEIGH_OUT"
  echo "Case directory: ${CASE_DIR}" >> "$NEIGH_OUT"
  echo "==================================================" >> "$NEIGH_OUT"

  echo "[${PX}x${PY}] Generating GridDist in ${CASE_DIR} ..." | tee -a "$LOG_OUT"
  (
    cd "$CASE_DIR"
    ../../GridDist $PX $PY $SIZE $SIZE
  ) | tee -a "$LOG_OUT"
  echo "" | tee -a "$LOG_OUT"

  echo "[${PX}x${PY}] Reading neighbour counts from input files..." | tee -a "$LOG_OUT"
  for r in $(seq 0 $((NP-1))); do
    f="${CASE_DIR}/input${NP}-${r}.dat"
    n=$(grep -m1 "^Neighbours:" "$f" | awk '{print $2}')
    printf "rank %2d: Neighbours = %s\n" "$r" "$n" >> "$NEIGH_OUT"
  done
  echo "" >> "$NEIGH_OUT"

  echo "[${PX}x${PY}] Optional solver run..." | tee -a "$LOG_OUT"
  if [ "$NP" -le 4 ]; then
    srun --nodes=1 --ntasks=$NP --ntasks-per-node=$NP --chdir="$CASE_DIR" ../../MPI_Fempois \
      | tee -a "$LOG_OUT"
  else
    srun --nodes=2 --ntasks=$NP --distribution=block --chdir="$CASE_DIR" ../../MPI_Fempois \
      | tee -a "$LOG_OUT"
  fi
  echo "" | tee -a "$LOG_OUT"
}

run_case 2 2
run_case 3 3

echo "Done." | tee -a "$LOG_OUT"
echo "Neighbour counts written to: $NEIGH_OUT" | tee -a "$LOG_OUT"
echo "Logs written to: $LOG_OUT" | tee -a "$LOG_OUT"

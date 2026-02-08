#!/bin/bash
#SBATCH --job-name=FEM_Ex4_2
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

set -euo pipefail

# ------------------------
# CONFIG
# ------------------------
PX=2
PY=2
NP=$((PX*PY))
SIZES=(100 200 400)

CASE_ROOT="data_ex_4_2"
mkdir -p "$CASE_ROOT"
OUT_TXT="${CASE_ROOT}/ex4_2_timebreakdown.txt"
OUT_CSV="${CASE_ROOT}/ex4_2_timebreakdown.csv"

# Start fresh
> "$OUT_TXT"
echo "grid,iters,t_calc,t_exch,t_glob,t_idle,ratio" > "$OUT_CSV"

echo "==================================================" | tee -a "$OUT_TXT"
echo "Exercise 4.2 FEM timing breakdown (P=${NP}, grid=${PX}x${PY})" | tee -a "$OUT_TXT"
echo "Output TXT: ${OUT_TXT}" | tee -a "$OUT_TXT"
echo "Output CSV: ${OUT_CSV}" | tee -a "$OUT_TXT"
echo "==================================================" | tee -a "$OUT_TXT"
echo "" | tee -a "$OUT_TXT"

# ------------------------
# BUILD
# ------------------------
echo "Compiling..." | tee -a "$OUT_TXT"
make | tee -a "$OUT_TXT"
echo "" | tee -a "$OUT_TXT"

# ------------------------
# RUNS
# ------------------------
for N in "${SIZES[@]}"; do
  CASE_DIR="${CASE_ROOT}/case_${PX}x${PY}_N${N}"
  mkdir -p "$CASE_DIR"

  echo "----------------------------------------" | tee -a "$OUT_TXT"
  echo "Case: ${PX}x${PY}, grid ${N}x${N}" | tee -a "$OUT_TXT"
  echo "Case directory: ${CASE_DIR}" | tee -a "$OUT_TXT"
  echo "----------------------------------------" | tee -a "$OUT_TXT"

  # MPI_Fempois expects input.dat (and sources.dat for GridDist) in CWD
  cp -f sources.dat "${CASE_DIR}/sources.dat"
  cp -f input.dat   "${CASE_DIR}/input.dat"

  echo "Generating GridDist..." | tee -a "$OUT_TXT"
  ( cd "$CASE_DIR" && ../../GridDist "$PX" "$PY" "$N" "$N" ) | tee -a "$OUT_TXT"
  echo "" | tee -a "$OUT_TXT"

  echo "Running MPI_Fempois..." | tee -a "$OUT_TXT"
  (
    cd "$CASE_DIR"
    srun --ntasks="$NP" ../../MPI_Fempois 2>&1 | tee run.log
  ) | tee -a "$OUT_TXT"
  echo "" | tee -a "$OUT_TXT"

  RUN_LOG="${CASE_DIR}/run.log"

  # Iterations: "Number of iterations : 141" -> field 5
  ITERS=$(grep -m1 -F "Number of iterations" "$RUN_LOG" | awk '{print $5}')

  # FIX: match the literal "(0) Time Breakdown:" using -F
  LINE0=$(grep -m1 -F "(0) Time Breakdown:" "$RUN_LOG" || true)
  if [[ -z "${LINE0}" ]]; then
    echo "WARNING: Could not find rank-0 Time Breakdown line for N=${N}" | tee -a "$OUT_TXT"
    continue
  fi

  # Extract numbers (strip trailing 's')
  t_calc=$(echo "$LINE0" | awk -F'Calc: ' '{print $2}' | awk -F's' '{print $1}')
  t_exch=$(echo "$LINE0" | awk -F'Exchange: ' '{print $2}' | awk -F's' '{print $1}')
  t_glob=$(echo "$LINE0" | awk -F'Global: ' '{print $2}' | awk -F's' '{print $1}')
  t_idle=$(echo "$LINE0" | awk -F'Idle: ' '{print $2}' | awk -F's' '{print $1}')

  ratio=$(awk -v c="$t_calc" -v e="$t_exch" -v g="$t_glob" 'BEGIN{printf "%.6f", c/(e+g)}')

  echo "${N},${ITERS},${t_calc},${t_exch},${t_glob},${t_idle},${ratio}" >> "$OUT_CSV"
done

echo "==================================================" | tee -a "$OUT_TXT"
echo "Done." | tee -a "$OUT_TXT"
echo "CSV written to: ${OUT_CSV}" | tee -a "$OUT_TXT"
echo "==================================================" | tee -a "$OUT_TXT"

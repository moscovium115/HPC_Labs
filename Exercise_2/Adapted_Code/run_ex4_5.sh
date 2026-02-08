#!/bin/bash
#SBATCH --job-name=FEM_Ex4_5
#SBATCH --partition=compute-p1
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2025
module load openmpi

set -euo pipefail

# ------------------------
# OUTPUT DIRS / FILES
# ------------------------
CASE_ROOT="data_ex_4_5"
mkdir -p "$CASE_ROOT"

OUT_TXT="${CASE_ROOT}/ex4_5_measurements.txt"
OUT_CSV="${CASE_ROOT}/ex4_5_measurements.csv"

# Start fresh
> "$OUT_TXT"
echo "part,px,py,p,n,iters,t_calc,t_exch,t_glob,t_idle,t_comm,t_comm_plus_idle,ratio,ratio_plus_idle" > "$OUT_CSV"

echo "==================================================" | tee -a "$OUT_TXT"
echo "Exercise 4.5 FEM measurements" | tee -a "$OUT_TXT"
echo "TXT: ${OUT_TXT}" | tee -a "$OUT_TXT"
echo "CSV: ${OUT_CSV}" | tee -a "$OUT_TXT"
echo "==================================================" | tee -a "$OUT_TXT"
echo "" | tee -a "$OUT_TXT"

# ------------------------
# BUILD
# ------------------------
echo "Compiling..." | tee -a "$OUT_TXT"
make | tee -a "$OUT_TXT"
echo "" | tee -a "$OUT_TXT"

# ------------------------
# HELPERS
# ------------------------
extract_rank0_breakdown() {
  local run_log="$1"
  grep -m1 -F "(0) Time Breakdown:" "$run_log" || true
}

extract_iters() {
  local run_log="$1"
  grep -m1 -F "Number of iterations" "$run_log" | awk '{print $5}' || true
}

run_case() {
  local part="$1"
  local px="$2"
  local py="$3"
  local n="$4"

  local p=$((px*py))
  local case_dir="${CASE_ROOT}/case_${part}_P${p}_${px}x${py}_N${n}"
  mkdir -p "$case_dir"

  echo "----------------------------------------" | tee -a "$OUT_TXT"
  echo "Case: part=${part}, P=${p} (${px}x${py}), grid=${n}x${n}" | tee -a "$OUT_TXT"
  echo "Dir:  ${case_dir}" | tee -a "$OUT_TXT"
  echo "----------------------------------------" | tee -a "$OUT_TXT"

  cp -f sources.dat "${case_dir}/sources.dat"
  cp -f input.dat   "${case_dir}/input.dat"

  echo "Generating GridDist..." | tee -a "$OUT_TXT"
  ( cd "$case_dir" && ../../GridDist "$px" "$py" "$n" "$n" ) | tee -a "$OUT_TXT"
  echo "" | tee -a "$OUT_TXT"

  echo "Running MPI_Fempois..." | tee -a "$OUT_TXT"
  local run_log="${case_dir}/run.log"

  # write solver output directly to run.log (reliable parsing)
  ( cd "$case_dir" && srun --ntasks="$p" ../../MPI_Fempois ) >"$run_log" 2>&1

  cat "$run_log" | tee -a "$OUT_TXT"
  echo "" | tee -a "$OUT_TXT"

  local line iters
  line="$(extract_rank0_breakdown "$run_log")"
  iters="$(extract_iters "$run_log")"

  if [[ -z "$line" ]]; then
    echo "WARNING: Missing rank-0 breakdown line in ${run_log}" | tee -a "$OUT_TXT"
    return 0
  fi
  if [[ -z "$iters" ]]; then
    echo "WARNING: Missing iteration count in ${run_log}" | tee -a "$OUT_TXT"
    return 0
  fi

  local t_calc t_exch t_glob t_idle
  t_calc=$(echo "$line" | awk -F'Calc: '     '{print $2}' | awk -F's' '{print $1}')
  t_exch=$(echo "$line" | awk -F'Exchange: ' '{print $2}' | awk -F's' '{print $1}')
  t_glob=$(echo "$line" | awk -F'Global: '   '{print $2}' | awk -F's' '{print $1}')
  t_idle=$(echo "$line" | awk -F'Idle: '     '{print $2}' | awk -F's' '{print $1}')

  local t_comm t_comm_plus_idle ratio ratio_plus_idle
  t_comm=$(awk "BEGIN{printf \"%.6f\", ${t_exch}+${t_glob}}")
  t_comm_plus_idle=$(awk "BEGIN{printf \"%.6f\", ${t_exch}+${t_glob}+${t_idle}}")
  ratio=$(awk "BEGIN{printf \"%.6f\", ${t_calc}/(${t_exch}+${t_glob})}")
  ratio_plus_idle=$(awk "BEGIN{printf \"%.6f\", ${t_calc}/(${t_exch}+${t_glob}+${t_idle})}")

  echo "${part},${px},${py},${p},${n},${iters},${t_calc},${t_exch},${t_glob},${t_idle},${t_comm},${t_comm_plus_idle},${ratio},${ratio_plus_idle}" >> "$OUT_CSV"
}

# ------------------------
# EXPERIMENTS
# ------------------------

# Part A: fixed P=4 (2x2)
PX_A=2
PY_A=2
SIZES_A=(50 100 150 200 300 400)

for n in "${SIZES_A[@]}"; do
  run_case "A_fixedP4" "$PX_A" "$PY_A" "$n"
done

# Part B: fixed n=1000, vary P (keep small under 2 minutes!)
N_B=1000
P_LIST_B=("1 1" "1 2")

for cfg in "${P_LIST_B[@]}"; do
  px=$(echo "$cfg" | awk '{print $1}')
  py=$(echo "$cfg" | awk '{print $2}')
  run_case "B_fixedN1000" "$px" "$py" "$N_B"
done

echo "==================================================" | tee -a "$OUT_TXT"
echo "Done." | tee -a "$OUT_TXT"
echo "CSV written to: ${OUT_CSV}" | tee -a "$OUT_TXT"
echo "==================================================" | tee -a "$OUT_TXT"

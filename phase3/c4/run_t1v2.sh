#!/usr/bin/env bash
# Run a stream of T1-v2 jobs. Usage: run_t1v2.sh <out_dir> "<model> <L> <level> <lres>" ...
# Records exact command, exit code and UTC start/end per run in <out_dir>/RUNS.txt.
set -u
OUT="$1"; shift
mkdir -p "$OUT/logs"
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
PY=${PY:-/home/claude/venv312/bin/python}
for spec in "$@"; do
  set -- $spec; MO=$1; L=$2; LEV=$3; LR=$4
  tag="${MO}_lres${LR}_L${L}_lev${LEV}"
  st=$(date -u +%FT%T.%3NZ)
  $PY phase3/c4/t1v2.py --model "$MO" --L "$L" --level "$LEV" --lres "$LR" --out "$OUT" > "$OUT/logs/t1v2_${tag}.log" 2>&1
  ec=$?
  echo "model=$MO L=$L level=$LEV lres=$LR exit=$ec start=$st end=$(date -u +%FT%T.%3NZ) cmd='python phase3/c4/t1v2.py --model $MO --L $L --level $LEV --lres $LR --out $OUT'" >> "$OUT/RUNS.txt"
done

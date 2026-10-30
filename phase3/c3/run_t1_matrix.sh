#!/usr/bin/env bash
# Run one stream of the T1 matrix. Usage: run_t1_matrix.sh <out_dir> "<L> <level> <l3>" ...
# Records exact command, exit code and UTC start/end per run in <out_dir>/RUNS.txt.
set -u
OUT="$1"; shift
mkdir -p "$OUT/logs"
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
PY=/home/claude/venv312/bin/python
for spec in "$@"; do
  set -- $spec; L=$1; LEV=$2; L3=$3
  st=$(date -u +%FT%T.%3NZ)
  $PY phase3/c3/t1.py --L "$L" --level "$LEV" --l3 "$L3" --out "$OUT" > "$OUT/logs/t1_L${L}_lev${LEV}.log" 2>&1
  ec=$?
  echo "L=$L level=$LEV l3=$L3 exit=$ec start=$st end=$(date -u +%FT%T.%3NZ) cmd='python phase3/c3/t1.py --L $L --level $LEV --l3 $L3 --out $OUT'" >> "$OUT/RUNS.txt"
done

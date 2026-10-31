#!/usr/bin/env bash
# Generic recorded runner. Usage: run_jobs.sh <out_dir> <script> "<args...>" ...
# Each argument string is passed to the script together with --out <out_dir>.
set -u
OUT="$1"; SCRIPT="$2"; shift 2
mkdir -p "$OUT/logs"
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
PY=${PY:-/home/claude/venv312/bin/python}
for args in "$@"; do
  tag=$(echo "$args" | tr -s ' -' '_' | sed 's/^_//')
  st=$(date -u +%FT%T.%3NZ)
  $PY "$SCRIPT" $args --out "$OUT" > "$OUT/logs/${tag}.log" 2>&1
  ec=$?
  echo "args='$args' exit=$ec start=$st end=$(date -u +%FT%T.%3NZ) cmd='python $SCRIPT $args --out $OUT'" >> "$OUT/RUNS.txt"
done

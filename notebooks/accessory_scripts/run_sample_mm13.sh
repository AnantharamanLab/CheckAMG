#!/bin/bash
# Per-sample runner of the DRAM-V mm13 re-run, usage: run_sample_mm13.sh SUBSET SAMPLE THREADS
set -euo pipefail
SUBSET="$1"; SAMPLE="$2"; THREADS="$3"
: "${DATA_ROOT:?set DATA_ROOT to the benchmark data root}"
: "${TMP_ROOT:?set TMP_ROOT to the tmp root}"
: "${DRAM_ENV:?set DRAM_ENV to the prefix of the DRAM_1.5.0_mm13 environment}"
# Directories appended to PATH after the DRAM env, the run used the peatlands_env and base conda bin directories for pigz
EXTRA_PATH="${EXTRA_PATH:-}"
WORK=$TMP_ROOT/dramv_mm13_rerun
IN=$DATA_ROOT/virsorter_outputs/$SUBSET/$SAMPLE/for-dramv
OUT=$DATA_ROOT/dramv_outputs_mm13/$SUBSET/$SAMPLE
ANN=$DATA_ROOT/dramv_annotations_mm13/$SUBSET/$SAMPLE
ARCH=$WORK/working_dir_archives/$SUBSET
LOGD=$WORK/logs/$SUBSET
export PATH=$DRAM_ENV/bin:$PATH${EXTRA_PATH:+:$EXTRA_PATH} TMPDIR=$TMP_ROOT OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore
command -v pigz >/dev/null || { echo "pigz not found" >&2; exit 9; }
mkdir -p "$LOGD" "$(dirname "$OUT")" "$ANN" "$ARCH"
[ "$(ls "$IN"/*.fa | wc -l)" -eq 1 ] && [ "$(ls "$IN"/*.tab | wc -l)" -eq 1 ]
FA=$(ls "$IN"/*.fa); TAB=$(ls "$IN"/*.tab)
if [ ! -e "$OUT/distilled/amg_summary.tsv" ]; then
  [ -e "$OUT" ] && mv "$OUT" "$OUT.incomplete.$(date +%Y%m%d%H%M%S)"
  { echo "start $(date -Is) threads=$THREADS"; mmseqs version; python -c 'import importlib.metadata as m; print("dram", m.version("DRAM-bio"))'; } > "$LOGD/${SAMPLE}_versions.txt"
  /usr/bin/time -v DRAM-v.py annotate -i "$FA" -v "$TAB" -o "$OUT" --keep_tmp_dir --threads "$THREADS" > "$LOGD/${SAMPLE}_annotate.log" 2> "$LOGD/${SAMPLE}_annotate.time"
  /usr/bin/time -v DRAM-v.py distill -i "$OUT/annotations.tsv" -o "$OUT/distilled" --max_auxiliary_score 4 --log_file_path "$OUT/distilled/distill.log" > "$LOGD/${SAMPLE}_distill.log" 2> "$LOGD/${SAMPLE}_distill.time"
  echo "dram_end $(date -Is)" >> "$LOGD/${SAMPLE}_versions.txt"
fi
if [ -d "$OUT/working_dir" ]; then
  TMPB="$OUT/working_dir/final-viral-combined-for-dramv/tmp"
  N=$(ls "$TMPB"/*.b6 | wc -l); [ "$N" -eq 8 ]
  for f in "$TMPB"/*.b6; do cp -n "$f" "$ANN/"; cmp "$f" "$ANN/$(basename "$f")"; done
  ( cd "$OUT" && tar -cf - working_dir | pigz -p "$THREADS" > "$ARCH/${SAMPLE}_working_dir.tar.gz.part" )
  pigz -t -p "$THREADS" "$ARCH/${SAMPLE}_working_dir.tar.gz.part"
  NT=$(pigz -dc -p "$THREADS" "$ARCH/${SAMPLE}_working_dir.tar.gz.part" | tar -tf - | grep -vc '/$' || true)
  NF=$(find "$OUT/working_dir" -type f | wc -l)
  [ "$NT" -eq "$NF" ]
  mv "$ARCH/${SAMPLE}_working_dir.tar.gz.part" "$ARCH/${SAMPLE}_working_dir.tar.gz"
  rm -rf "$OUT/working_dir"
  echo "archived $(date -Is) files=$NF b6=$N" >> "$LOGD/${SAMPLE}_versions.txt"
fi
echo "done $(date -Is)" >> "$LOGD/${SAMPLE}_versions.txt"

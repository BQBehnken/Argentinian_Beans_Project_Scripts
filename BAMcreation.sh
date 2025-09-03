#!/usr/bin/env bash
set -euo pipefail

# ---- file path ----------------------------------------------------
REF="/media/.../Pvulgaris_442_v2.0.fa"
RAWDIR="/media/.../raw_data"
THREADS="${THREADS:-40}" # because we have 40 
# -----------------------------------------------------------------------------

# find all *_1.fq (uncompressed)
mapfile -t R1S < <(find "$RAWDIR" -type f -name "*_1.fq" | sort)

if (( ${#R1S[@]} == 0 )); then
  echo "[ERROR] No *_1.fq files found under $RAWDIR" >&2
  exit 1
fi

for r1 in "${R1S[@]}"; do
  r2="${r1/_1.fq/_2.fq}"
  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing mate for: $r1 (expected $r2). Skipping." >&2
    continue
  fi

  # Sample naming for internal alias: D266_* -> 266_v2
  bn="$(basename "$r1")"                # e.g., D266_USPD..._L7_1.fq
  sample_tag="${bn%%_*}"                # D266
  numeric="${sample_tag#D}"             # 266
  out="${numeric}_v2"                   # 266_v2 # we recreated the BAMs as a sanity check, so the v2 convention stuck

  echo "[INFO] Mapping $sample_tag -> $out"

  bwa mem "$REF" "$r1" "$r2" -o "${out}.sam" -t "$THREADS"
  samtools view -@ "$THREADS" -bS "${out}.sam" > "${out}.bam"
  samtools sort -@ "$THREADS" "${out}.bam" -o "${out}.sorted.bam"
  samtools index -@ "$THREADS" "${out}.sorted.bam"
done

echo "[DONE] Mapped, sorted, and indexed all FASTQ pairs found."

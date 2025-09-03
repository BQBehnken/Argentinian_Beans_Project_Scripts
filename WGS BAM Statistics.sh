#!/usr/bin/env bash
set -euo pipefail

# --- sanity checks ---
for tool in mosdepth samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[!] $tool not found in PATH"; exit 1; }
done

# BAMs
pv_samples=(Accesion1 Accession2 ... Accession21) # all my accession numbers in folder

run_for_bam () {
  local prefix="$1"   # mosdepth output prefix
  local bam="$2"      # path to BAM file

  if [[ -f "$bam" ]]; then
    # ensure index exists for mosdepth
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
      echo "[*] Index missing for $bam — indexing..."
      samtools index -@ 8 "$bam"
    fi

    echo "[*] Running mosdepth on $bam"
    mosdepth "$prefix" "$bam"

    echo "[*] Running samtools stats for $bam ..."
    samtools stats -@ 8 "$bam" > "${bam%.bam}.samstats.txt"
  else
    echo "[!] BAM file not found: $bam" >&2
  fi
}

# BAM is sample.bam
for s in "${pv_samples[@]}"; do
  run_for_bam "$s" "${s}.bam"
done

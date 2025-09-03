#!/usr/bin/env bash
set -euo pipefail

# ---- Edit these to your setup ------------------------------------------------
REF="Pvulgaris_442_v2.0.fa"             # full reference FASTA (in current dir)
REGION="Chr07:7407744-7413510"          # INR locus (must match BAM naming)
BAM_DIR="."                              # directory with BAMs ('.' = current)
GLOB="*.bam"                             # which BAMs to process
OUT_VCF="results/vcf"                    # where to write VCFs
OUT_FA="results/consensus"               # where to write consensus FASTAs
LOGS="logs"                              # per-sample logs
# Haplotype policy for consensus on bcftools:
#   -H 1 or -H 2 (choose hap1/hap2), or -I for IUPAC codes.
CONSENSUS_ARGS="-H 1"
# -----------------------------------------------------------------------------

mkdir -p "$OUT_VCF" "$OUT_FA" "$LOGS"

# 0) Ensure the full REF has a .fai index
if [[ ! -f "${REF}.fai" ]]; then
  echo "[INFO] Indexing REF with samtools faidx"
  samtools faidx "$REF"
fi

# 1) Make the trimmed INR reference once
INR_REF="INR_reference442.fa"
echo "[INFO] Creating trimmed reference $INR_REF for $REGION"
samtools faidx "$REF" "$REGION" > "$INR_REF"

# quick sanity check: header should print
head -1 "$INR_REF" | grep -q ">" || { echo "[ERROR] Trimmed reference looks wrong"; exit 1; }

# 2) Process each BAM
shopt -s nullglob
for bam in "$BAM_DIR"/$GLOB; do
  sample="$(basename "$bam" .bam)"
  echo "[INFO] Processing $sample"

  # 2a) Ensure BAM is indexed
  [[ -f "${bam}.bai" ]] || samtools index "$bam"

  # 2b) mpileup -> call -> bgzip VCF  (restricted to INR REGION)
  bcftools mpileup -Ou -f "$REF" -r "$REGION" "$bam" 2> "$LOGS/${sample}.log" \
  | bcftools call   -mv -Oz -o "$OUT_VCF/${sample}.vcf.gz" 2>> "$LOGS/${sample}.log"

  # 2c) index the VCF
  tabix -p vcf "$OUT_VCF/${sample}.vcf.gz"

  # 2d) consensus against the **trimmed** INR reference (no -r here)
  #     (choose haplotype with -H 1 / -H 2 or use -I for IUPAC on your version)
  bcftools consensus -f "$INR_REF" $CONSENSUS_ARGS \
    "$OUT_VCF/${sample}.vcf.gz" > "$OUT_FA/${sample}.fa"
   
done

  # 3a) Merge all per-sample VCFs into one multi-sample VCF
  
if compgen -G "$OUT_VCF/*.vcf.gz" > /dev/null; then
  echo "[INFO] Merging $(ls "$OUT_VCF"/*.vcf.gz | wc -l) VCF(s)"
  bcftools merge "$OUT_VCF"/*.vcf.gz -Oz -o "$OUT_VCF/inr_merged.vcf.gz"
  tabix -p vcf "$OUT_VCF/inr_merged.vcf.gz"
  echo "[INFO] Merged VCF written to: $OUT_VCF/inr_merged.vcf.gz"
else
  echo "[WARN] No VCFs found to merge in $OUT_VCF"
fi

  # 3b) Build a single multi-FASTA from all per-sample consensus FASTAs
if compgen -G "$OUT_FA/*.fa" > /dev/null; then
  cat $(ls "$OUT_FA"/*.fa | sort) > "$OUT_FA/inr_consensus.multi.fa"
  echo "[INFO] Multi-FASTA: $OUT_FA/inr_consensus.multi.fa"
else
  echo "[WARN] No FASTAs found to merge in $OUT_FA"
fi
echo "[DONE] Consensus FASTAs in $OUT_FA ; VCFs in $OUT_VCF"

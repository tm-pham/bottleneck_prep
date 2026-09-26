#!/usr/bin/env bash
# Build a bottleneck_prep samplesheet from an mtb-call2 results directory.
#
#   bin/make_samplesheet.sh /path/to/mtb-call2/results > samplesheet.tsv
#
# Expects mtb-call2's layout: results/<batch>/<sample>/{bams,vars}/
# and that BAM publishing was switched on in mtb-call2 (see README.md).
set -euo pipefail
results=${1:?usage: make_samplesheet.sh <mtb-call2 results dir>}

printf 'sample_id\tbam\tlofreq_vcf\n'
for dir in "$results"/*/*/; do
  sample=$(basename "$dir")
  bam=$(ls "$dir"bams/"${sample}"_filtered.bam 2>/dev/null || ls "$dir"bams/"${sample}"_merged_mrkdup.bam 2>/dev/null || true)
  vcf=$(ls "$dir"vars/"${sample}"_lofreq_unfilt.vcf.gz 2>/dev/null || true)
  if [[ -z "$bam" || -z "$vcf" ]]; then
    echo "WARNING: skipping $sample (missing BAM or unfiltered LoFreq VCF)" >&2
    continue
  fi
  printf '%s\t%s\t%s\n' "$sample" "$(realpath "$bam")" "$(realpath "$vcf")"
done

# bottleneck_prep

A small Nextflow pipeline that turns **mtb-call2 output** (per-sample BAM + unfiltered LoFreq VCF) and a **hand-verified donor–recipient pair list** into one table of allele counts per pair × site. That table is the input to the beta-binomial estimator in `code/simulation/R/estimate.R`.

This is a first pass, kept simple on purpose. Per-sample processing (trimming, Kraken2, mapping, LoFreq) stays in [ksw9/mtb-call2](https://github.com/ksw9/mtb-call2). This pipeline only does the per-pair step that mtb-call2 doesn't cover.

```
FASTQ ──mtb-call2──► BAM + LoFreq VCF (per sample)
                                │
            pairs.tsv ──────────┤
                                ▼
                      bottleneck_prep (this)
          COUNT_PAIR (per pair)  →  BUILD_BB_TABLE
                                ▼
            bb_sites.tsv, bb_pair_summary.tsv
                                ▼
          R/as_bb_input.R  →  estimate.R::fit_bb()
```

## What it does

For each pair:

1. **Candidate sites.** Takes the union of SNPs that LoFreq called (unfiltered) in the donor *or* the recipient.
2. **Re-counts alleles at every candidate site in both BAMs** with `bcftools mpileup` (MAPQ ≥ 20, BQ ≥ 20, no BAQ, no depth cap). A donor iSNV missing from the recipient is recorded as `k = 0` at the recipient's real depth, not dropped. For the beta-binomial, those lost variants are the informative sites.
3. **Flags, doesn't drop,** sites in the PE/PPE mask (`assets/H37Rv_ppe.bed.gz`, the same BED as mtb-call2).

No MAF threshold, depth filter or mask is applied in the pipeline. They are arguments of `as_bb_input()` in R, so you can sweep them (e.g. T = 0.5/1/2/5%) without re-running Nextflow.

## Step 0: settings to change in mtb-call2

mtb-call2 doesn't publish BAMs by default. In `mtb-call2/main.nf`, uncomment the two `bam_files` lines: one in `publish:` and one in the `output {}` block. Keep `run_lofreq = true` (the default). LoFreq runs per sample there, which avoids the joint-calling circularity described in `literature/TB_bottleneck_pipeline_design.md` §6.4.

## Inputs

**samplesheet.tsv**, one row per sample. `bai` is optional and defaults to `<bam>.bai`.

| sample_id | bam | lofreq_vcf |
|---|---|---|
| S1 | /path/S1_filtered.bam | /path/S1_lofreq_unfilt.vcf.gz |

To generate it from an mtb-call2 results folder: `bin/make_samplesheet.sh /path/to/mtb-call2/results > samplesheet.tsv`

**pairs.tsv**, one row per pair, curated by hand *before* looking at sequence data. It needs `pair_id`, `donor_sample_id` and `recipient_sample_id`; any extra columns (cohort, dates, link type, …) are carried into the output. See `assets/pairs_template.tsv`.

**reference**: the H37Rv FASTA used by mtb-call2, with its `.fai`.

## Running

```bash
# on CHPC (SLURM + Apptainer; reuses the mtb-call2 variantcalling image)
module load nextflow apptainer
nextflow run main.nf -profile chpc \
  --samplesheet samplesheet.tsv \
  --pairs pairs.tsv \
  --reference /path/to/mtb-call2/resources/refs/H37Rv.fasta \
  --outdir results/colangeli

# built-in synthetic test (needs bwa, samtools, bgzip, bcftools, Rscript on PATH)
python3 tests/make_test_data.py
nextflow run main.nf -profile test
```

## Outputs

- `bb_sites.tsv`: one row per pair × site × alt allele. It contains:
  - the pair metadata;
  - `called_in_donor` / `called_in_recipient`;
  - LoFreq AF/QUAL/SB for each host;
  - re-counted `*_ref`, `*_alt`, `*_depth` and `*_af` for each host;
  - strand-split alt counts (`*_alt_fwd`, `*_alt_rev`);
  - `mask_pe_ppe`.
- `bb_pair_summary.tsv`: sites per pair (total, donor-called, recipient-called, shared, donor-called outside the mask) and median depths. Check this first: if most pairs have almost no unmasked donor iSNVs, nothing downstream will fix that.
- `pairs/<pair_id>/`: the raw per-pair LoFreq extracts and mpileup counts, kept for debugging.

## Into the estimator

```r
source("code/simulation/R/estimate.R")
source("code/bottleneck_prep/R/as_bb_input.R")

sites <- read_bb_sites("code/bottleneck_prep/results/colangeli/bb_sites.tsv")
dat   <- as_bb_input(sites, T = 0.01, min_depth = 20,
                     min_donor_af = 0.01, max_donor_af = 0.99,
                     exclude_pe_ppe = TRUE)
bb_site_counts(dat)          # sites per pair after filtering
fit <- fit_bb(dat, T = 0.01) # use the SAME T in both calls
summarise_fit(fit)
```

`as_bb_input()` does the following:

- keeps donor-called sites only (the Sobel Leonard definition);
- requires `min_depth` in **both** hosts, dropping failing sites rather than recording them as zero;
- sets recipient counts below T to 0 so `fit_bb()`'s threshold correction applies;
- uses the re-counted donor frequency as `nu_D`.

## Known simplifications (first pass)

- SNPs only. Indels called by LoFreq are ignored.
- The PE/PPE mask is the only mask. No repeat, IS-element or phage masks yet.
- Frequencies are recounted without BAQ, while LoFreq uses BAQ, so `donor_af` can differ slightly from `donor_lofreq_af`. Both are kept.
- No mixed-infection or contamination screen at the pair level. Use mtb-call2's Kraken2 and TB-Profiler outputs to exclude samples before building `pairs.tsv`.
- Output is TSV, not Parquet.

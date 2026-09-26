# *M. tuberculosis* transmission bottleneck data preparation

Pipeline to prepare within-host variant data from epidemiologically linked donor–recipient pairs for transmission bottleneck estimation. It takes per-sample BAMs and unfiltered LoFreq VCFs produced by [mtb-call2](https://github.com/ksw9/mtb-call2), re-counts alleles at every candidate site in both members of each pair, and outputs a single table ready for beta-binomial bottleneck estimators. The pipeline is written in Nextflow.

## Pipeline steps

* **COUNT_PAIR** (per pair)

  * Collects SNPs called by LoFreq in the donor or the recipient.

  * Re-counts reference/alternate (and A/C/G/T) reads at every site in **both** samples with bcftools mpileup, so variants lost in transmission are kept as zero counts at the real depth.

* **BUILD_BB_TABLE** (all pairs)

  * Combines pairs into one table and flags sites in PE/PPE regions (flagged, not removed).

  * Writes a per-pair summary of site counts and depth.

No frequency threshold, depth filter or mask is applied in the pipeline; these are set in the R step so they can be varied without re-running Nextflow.

## Installation & set-up

1. Run [mtb-call2](https://github.com/ksw9/mtb-call2) on your samples with `run_lofreq = true`.\
N.B. mtb-call2 does not keep BAMs by default: uncomment the `bam_files` lines in both the `publish:` and `output {}` sections of its **main.nf**.

2. Clone this repository.

```
git clone https://github.com/tm-pham/bottleneck_prep.git
```

3. Load Nextflow and your container tool (Apptainer/Singularity or Docker).

```
module load nextflow apptainer
```

## Running analyses

1. Create a tab-delimited **samplesheet**. It can be generated from an mtb-call2 results directory:

```
bin/make_samplesheet.sh /path/to/mtb-call2/results > samplesheet.tsv
```

| sample_id | bam | lofreq_vcf |
| :---: | :---: | :---: |
| unique_sample_identifier | /full/path/sample_filtered.bam | /full/path/sample_lofreq_unfilt.vcf.gz |

2. Create a tab-delimited **pairs** file, curated from the epidemiological data before looking at the sequences. Extra columns (cohort, dates, link type) are carried through to the output. See **assets/pairs_template.tsv**.

| pair_id | donor_sample_id | recipient_sample_id |
| :---: | :---: | :---: |
| unique_pair_identifier | donor_sample_identifier | recipient_sample_identifier |

3. Run the pipeline.

```
nextflow run main.nf -profile chpc \
  --samplesheet samplesheet.tsv \
  --pairs pairs.tsv \
  --reference /path/to/mtb-call2/resources/refs/H37Rv.fasta \
  --outdir results
```

N.B. Use `-profile docker` or `-profile local` outside CHPC. A synthetic test run is available with `python3 tests/make_test_data.py && nextflow run main.nf -profile test`.

4. Estimate the bottleneck size in R, e.g. with [ViralBottleneck](https://github.com/BowenArchaman/ViralBottleneck):

```r
source("R/as_bb_input.R")
source("R/viralbottleneck.R")
sites <- read_bb_sites("results/bb_sites.tsv")
run_viralbottleneck(sites, method = "Beta_binomial_Exact", variant_calling = 0.01, pooled = TRUE)
```

## Outputs

<pre>
<b>results</b>
│
├── <b>bb_sites.tsv</b>
│   One row per pair × site: read counts in donor and recipient, LoFreq annotations, PE/PPE flag
│
├── <b>bb_pair_summary.tsv</b>
│   Number of sites and median depth per pair
│
└── <b>pairs</b>
    └── <b>pair_id</b>
        Per-pair LoFreq extracts and allele counts
</pre>

See [docs/details.md](docs/details.md) for column definitions, estimator options and known limitations.

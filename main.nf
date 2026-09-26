#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * bottleneck_prep
 * ---------------
 * Turns per-sample mtb-call2 output (BAM + unfiltered LoFreq VCF) plus a
 * hand-curated donor-recipient pair list into one long table of allele counts
 * per pair x site, ready for the beta-binomial estimator (code/simulation/R/estimate.R).
 *
 * For every pair:
 *   1. candidate sites = union of SNPs LoFreq called in the donor OR the recipient
 *   2. ref/alt read counts are re-counted at every candidate site in BOTH samples
 *      (so a donor iSNV that is absent in the recipient is recorded as k = 0 at
 *      depth n, not silently missing)
 * Then all pairs are combined into bb_sites.tsv (+ a per-pair summary).
 *
 * No frequency threshold, depth filter or mask is applied here: those are
 * parameters of the R step (R/as_bb_input.R) so they can be swept without
 * re-running the pipeline.
 */

include { COUNT_PAIR     } from './modules/count_pair.nf'
include { BUILD_BB_TABLE } from './modules/build_bb_table.nf'

workflow {

    if (!params.samplesheet || !params.pairs || !params.reference) {
        error "Please provide --samplesheet, --pairs and --reference (see README.md)"
    }

    reference     = file(params.reference, checkIfExists: true)
    reference_fai = file("${params.reference}.fai", checkIfExists: true)
    mask_bed      = file(params.mask_bed, checkIfExists: true)
    pairs_file    = file(params.pairs, checkIfExists: true)

    // sample_id -> [bam, bai, vcf]
    samples = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { r ->
            def bai = r.bai ?: "${r.bam}.bai"
            tuple(r.sample_id, file(r.bam, checkIfExists: true),
                  file(bai, checkIfExists: true), file(r.lofreq_vcf, checkIfExists: true))
        }

    // pair_id, donor_id, recipient_id
    pairs = channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.donor_sample_id, r.pair_id, r.recipient_sample_id) }

    // attach donor files, then recipient files
    pair_inputs = pairs
        .combine(samples, by: 0)                                        // donor_id, pair, recip, d_bam, d_bai, d_vcf
        .map { d, p, r, dbam, dbai, dvcf -> tuple(r, p, d, dbam, dbai, dvcf) }
        .combine(samples, by: 0)                                        // recip, pair, donor, d..., r_bam, r_bai, r_vcf
        .map { r, p, d, dbam, dbai, dvcf, rbam, rbai, rvcf ->
            tuple(p, d, r, dbam, dbai, dvcf, rbam, rbai, rvcf) }

    COUNT_PAIR(pair_inputs, reference, reference_fai)

    BUILD_BB_TABLE(
        COUNT_PAIR.out.tables.map { it[1..3] }.flatten().collect(),
        pairs_file,
        mask_bed
    )
}

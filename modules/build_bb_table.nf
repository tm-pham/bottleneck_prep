process BUILD_BB_TABLE {

    // Combine all pairs into one long table (one row per pair x site x alt allele)
    // plus a per-pair summary. Base R only.

    label 'rscript'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tables
    path pairs
    path mask_bed

    output:
    path "bb_sites.tsv",        emit: sites
    path "bb_pair_summary.tsv", emit: summary

    script:
    """
    build_bb_table.R \\
        --pairs ${pairs} \\
        --mask ${mask_bed} \\
        --out-sites bb_sites.tsv \\
        --out-summary bb_pair_summary.tsv
    """

    stub:
    """
    touch bb_sites.tsv bb_pair_summary.tsv
    """
}

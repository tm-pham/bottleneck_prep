process COUNT_PAIR {

    // Per pair: union of LoFreq SNP sites in donor + recipient, then re-count
    // alleles at those sites in both BAMs with bcftools mpileup.

    tag "${pair_id}"
    label 'bcftools'
    publishDir "${params.outdir}/pairs/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), val(donor_id), val(recipient_id),
          path(donor_bam), path(donor_bai), path(donor_vcf),
          path(recipient_bam), path(recipient_bai), path(recipient_vcf)
    path reference
    path reference_fai

    output:
    tuple val(pair_id),
          path("${pair_id}.donor_lofreq.tsv"),
          path("${pair_id}.recipient_lofreq.tsv"),
          path("${pair_id}.counts.tsv"), emit: tables

    script:
    """
    # 1. LoFreq SNP calls (unfiltered) -> pair_id, chrom, pos, ref, alt, qual, af, dp, sb
    for role in donor recipient; do
        vcf=\$( [ "\$role" = donor ] && echo ${donor_vcf} || echo ${recipient_vcf} )
        bcftools view -v snps "\$vcf" \\
          | bcftools query -f '${pair_id}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AF\\t%INFO/DP\\t%INFO/SB\\n' \\
          > ${pair_id}.\${role}_lofreq.tsv
    done

    # 2. Union site set for this pair (CHROM, POS)
    cat ${pair_id}.donor_lofreq.tsv ${pair_id}.recipient_lofreq.tsv \\
      | cut -f2,3 | sort -k1,1 -k2,2n -u > sites.tsv

    # 3. Allele counts at every union site in BOTH samples.
    #    Column order is fixed by BAM order: donor first, recipient second.
    #    -B: no BAQ (raw base qualities, as for iSNV frequency estimation)
    #    -I: skip indels; -d: no practical depth cap
    printf 'pair_id\\tchrom\\tpos\\tref\\talts\\tdonor_AD\\tdonor_ADF\\tdonor_ADR\\trecipient_AD\\trecipient_ADF\\trecipient_ADR\\n' \\
      > ${pair_id}.counts.tsv
    if [ -s sites.tsv ]; then
        bcftools mpileup \\
            -f ${reference} -R sites.tsv \\
            --ignore-RG -B -I \\
            -q ${params.min_mapq} -Q ${params.min_baseq} \\
            -d ${params.max_depth} \\
            -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR \\
            ${donor_bam} ${recipient_bam} \\
          | bcftools query -f '${pair_id}\\t%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%ADF\\t%ADR]\\n' \\
          >> ${pair_id}.counts.tsv
    fi
    """

    stub:
    """
    touch ${pair_id}.donor_lofreq.tsv ${pair_id}.recipient_lofreq.tsv ${pair_id}.counts.tsv
    """
}

/*
 * Summarizes outputs for a particular sample-reference pair in .tsv format. 
 */
process BUILD_ALIGNMENT_SUMMARY {

    container 'community.wave.seqera.io/library/samtools_bc:9dec48711c1b0246'
    label 'process_low'

    input:
    tuple val(meta), path(merged_bam), path(consensus_bam), path(ref), val(ref_info), path(consensus_fa)

    output:
    tuple val(meta), path("*.tsv"),     emit: alignment_summary

    script:
    """
    # 1. sample info
    sample_name="${meta.id}"

    # 2. reference info
    ref_acc="${ref_info.acc}"
    ref_tag="${ref_info.tag}"
    ref_header="${ref_info.header}"

    # 3. mapping stats
    raw_reads=\$(samtools view -c ${merged_bam})
    mapped_reads=\$(samtools view -c -F 4 ${merged_bam})
    
    x100_reads_mapped=\$(echo "\$mapped_reads * 100" | bc)
    pct_reads_mapped=\$(echo "scale=2; \$x100_reads_mapped / \$raw_reads" | bc)

    # 4. coverage stats
    mean_depth=\$(samtools coverage ${merged_bam} | awk 'NR==2 {print \$7}')
    coverage=\$(samtools coverage ${merged_bam} | awk 'NR==2 {print \$6}')

    # build output .tsv
    {
        printf "sample_name\tref_acc\tref_tag\tref_header\traw_reads\tmapped_reads\tpct_reads_mapped\tmean_depth\tcoverage\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "\$sample_name" \
            "\$ref_acc" "\$ref_tag" "\$ref_header" \
            "\$raw_reads" "\$mapped_reads" "\$pct_reads_mapped" \
            "\$mean_depth" "\$coverage"
    } > "${ref.baseName}.tsv"
    """
}
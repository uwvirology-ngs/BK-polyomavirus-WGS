/*
 * Summarizes outputs for a particular sample-reference pair in .tsv format. 
 */
process BUILD_ALIGNMENT_SUMMARY {

    container 'community.wave.seqera.io/library/samtools_bc:9dec48711c1b0246'
    label 'process_low'

    input:
    tuple val(meta), path(merged_bam), path(consensus_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*.tsv"),     emit: alignment_summary

    script:
    """
    # 1. sample name
    sample_name="${meta.id}"

    # 2. reference info
    ref_acc="${ref_info.acc}"
    ref_tag="${ref_info.tag}"
    ref_header="${ref_info.header}"

    # 3. raw reads
    raw_reads=\$(samtools view -c ${merged_bam})
    
    # 4. mapping stats
    mapped_reads=\$(samtools view -c -F 4 ${merged_bam})
    x100_reads_mapped=\$(echo "\$mapped_reads * 100" | bc)
    pct_reads_mapped=\$(echo "scale=2; \$x100_reads_mapped / \$raw_reads" | bc)

    # 5. mean depth
    mean_depth=\$(samtools depth ${merged_bam} | awk '{accum+=\$3} END {print accum/NR}')

    # build output .tsv
    {
        printf "sample_name\tref_acc\tref_tag\tref_header\traw_reads\tmapped_reads\tpct_reads_mapped\tmean_depth\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "\$sample_name" \
            "\$ref_acc" "\$ref_tag" "\$ref_header" \
            "\$raw_reads" "\$mapped_reads" "\$pct_reads_mapped" \
            "\$mean_depth"
    } > "${ref.baseName}.tsv"
    """
}
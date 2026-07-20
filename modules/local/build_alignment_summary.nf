/*
 * Summarizes outputs for a particular sample-reference pair in .tsv format. 
 */
process BUILD_ALIGNMENT_SUMMARY {

    label 'twist'
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

    # 4. mean depth
    mean_depth=\$(samtools depth ${merged_bam} | awk '{accum+=\$3} END {print accum/NR}')

    # build output .tsv
    {
        printf "sample_name\tref_acc\tref_tag\tref_header\traw_reads\tmean_depth\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "\$sample_name" \
            "\$ref_acc" "\$ref_tag" "\$ref_header" \
            "\$raw_reads" \
            "\$mean_depth"
    } > "${ref.baseName}.tsv"
    """
}
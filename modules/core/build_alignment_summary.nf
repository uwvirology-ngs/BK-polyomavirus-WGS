/*
 * Summarizes outputs for a particular sample-reference pair in .tsv format. 
 */
process BUILD_ALIGNMENT_SUMMARY {

    container 'community.wave.seqera.io/library/samtools_bc:9dec48711c1b0246'
    label 'process_low'

    input:
    tuple val(meta), path(merged_bam), path(consensus_bam), path(ref), val(ref_info), path(consensus_fa), val(genomic_region), val(genomic_region_len)

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

    # 5. consensus genome stats
    consensus_bases=\$(grep -v '^>' ${consensus_fa} | tr -d '\\n')
    consensus_len=\${#consensus_bases}
    n_bases=\$(grep -o "N" <<< \$consensus_bases | wc -l)
    pct_ns=\$(echo "\${n_bases}/\${consensus_len}*100" | bc -l | awk 'FNR==1{print val,\$1}')
    pct_ns_format=\$(printf "%.4f" "\${pct_ns}")

    # 6. cds coverage
    samtools index ${merged_bam}
    pct_cds_covered=\$(samtools coverage -r ${genomic_region} ${merged_bam} | awk 'NR>1' | cut -f6)
    pct_cds_covered_formatted=\$(printf "%.2f" "\${pct_cds_covered}")
    mean_cds_coverage=\$(samtools coverage -r ${genomic_region} ${merged_bam} | awk 'NR>1' | cut -f7)
    mean_cds_coverage_formatted=\$(printf "%.2f" "\${mean_cds_coverage}")
    num_bases_cds_100x=\$(samtools depth -r ${genomic_region} ${merged_bam} | awk '{if(\$3>100)print\$3}' | wc -l)
    pct_cds_100x=\$(echo "\${num_bases_cds_100x}/${genomic_region_len}*100" | bc -l)
    pct_cds_100x_formatted=\$(printf "%.2f" "\${pct_cds_100x}")

    # build output .tsv
    {
        printf "sample_name\tref_acc\tref_tag\tref_header\traw_reads\tmapped_reads\tpct_reads_mapped\tmean_depth\tcoverage\tconsensus_len\tn_bases\tpct_ns\tpct_cds_covered\tmean_cds_coverage\tnum_bases_cds_100x\tpct_cds_100x\n"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "\$sample_name" \
            "\$ref_acc" "\$ref_tag" "\$ref_header" \
            "\$raw_reads" "\$mapped_reads" "\$pct_reads_mapped" \
            "\$mean_depth" "\$coverage" \
            "\$consensus_len" "\$n_bases" "\$pct_ns_format" \
            "\$pct_cds_covered_formatted" "\$mean_cds_coverage_formatted" "\$num_bases_cds_100x" "\$pct_cds_100x_formatted"
    } > "${ref.baseName}.tsv"
    """
}
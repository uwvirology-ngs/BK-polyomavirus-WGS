/*
 * An extra step for the Twist Bioscience UMI Protocol
 *
 * Builds a suite of metrics for quality control of duplex sequencing data
 * directly from GroupReadsByUMI output. 
 */
process FGBIO_COLLECT_DUPLEX_SEQ_METRICS {

    label 'twist'
    label 'process_medium'

    input:
    tuple val(meta), val(ref_info), path(ref), path(grouped_bam)

    output:
    tuple val(meta), val(ref_info), path("*"),  emit: duplex_seq_metrics

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    fgbio CollectDuplexSeqMetrics \\
        -Xmx${avail_mem}g \\
        --input=${grouped_bam} \\
        --output="${ref.baseName}" \\
        --duplex-umi-counts=true
    """
}

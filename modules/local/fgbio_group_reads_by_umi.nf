/*
 * Step 7 of the Twist Bioscience UMI Protocol
 *
 * Group reads which appear to map to the same molecule. 
 */
process FGBIO_GROUP_READS_BY_UMI {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(merged_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*.bam"), path(ref), val(ref_info),  emit: grouped_bam
    tuple val(meta), path("*histogram.txt"),                    emit: grouped_hist

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    fgbio GroupReadsByUmi \\
        -Xmx${avail_mem}g \\
        --strategy="paired" \\
        --input=${merged_bam} \\
        --output="${ref.simpleName}_grouped.bam" \\
        --raw-tag=RX \\
        --min-map-q=10 \\
        --edits=1 \\
        --family-size-histogram="${ref.simpleName}_histogram.txt"
    """
}
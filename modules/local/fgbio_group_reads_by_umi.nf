/*
 * Step 7 of the Twist Bioscience UMI Protocol
 *
 * Group reads which appear to map to the same molecule. 
 */
process FGBIO_GROUP_READS_BY_UMI {

    container 'community.wave.seqera.io/library/fgbio:1.5.1--9658b275f6ef481d'

    input:
    tuple val(meta), path(merged_bam), path(ref)

    output:
    tuple val(meta), path("*.bam"), path(ref),  emit: grouped_bam

    script:
    """
    fgbio GroupReadsByUmi \\
        --strategy="paired" \\
        --input=${merged_bam} \\
        --output="${ref.simpleName}_grouped.bam" \\
        --raw-tag=RX \\
        --min-map-q=10 \\
        --edits=1
    """
}
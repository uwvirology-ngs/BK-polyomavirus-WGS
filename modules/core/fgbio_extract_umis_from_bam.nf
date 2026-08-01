/*
 * Step 2 of the Twist Bioscience UMI Protocol
 *
 * Extract UMI sequences from the reads and store them as 
 * metadata within the unaligned BAM files. 
 */
process FGBIO_EXTRACT_UMIS_FROM_BAM {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(unaligned_bam)

    output:
    tuple val(meta), path("${meta.id}_umi_extracted.bam"),  emit: umi_extracted_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    fgbio ExtractUmisFromBam \\
        -Xmx${avail_mem}g \\
        --input="${unaligned_bam}" \\
        --output="${meta.id}_umi_extracted.bam" \\
        --read-structure=5M2S+T 5M2S+T \\
        --molecular-index-tags=ZA ZB \\
        --single-tag=RX
    """
}
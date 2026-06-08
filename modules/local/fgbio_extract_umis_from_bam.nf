/*
 * Step 2 of the Twist Bioscience UMI Protocol
 *
 * Extract UMI sequences from the reads and store them as 
 * metadata within the unaligned BAM files. 
 */
process FGBIO_EXTRACT_UMIS_FROM_BAM {

    container 'community.wave.seqera.io/library/fgbio:1.5.1--9658b275f6ef481d'

    input:
    tuple val(meta), path(unaligned_bam)

    output:
    tuple val(meta), path("${meta.id}_umi_extracted.bam"),  emit: umi_extracted_bam

    script:
    """
    fgbio ExtractUmisFromBam \\
        --input="${unaligned_bam}" \\
        --output="${meta.id}_umi_extracted.bam" \\
        --read-structure=5M2S+T 5M2S+T \\
        --molecular-index-tags=ZA ZB \\
        --single-tag=RX
    """
}
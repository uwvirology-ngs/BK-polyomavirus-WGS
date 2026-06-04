/*
 * With the UMI bases in the unaligned BAM file, convert the unaligned BAM to FastQ 
 * so we may select a reference then align the reads.
 * 
 * Does not use interleaved reads as it's not convenient for revica-strm.
 */
process PICARD_SAM_TO_FASTQ {

    container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

    input:
    tuple val(meta), path(unaligned_bam_umi_extracted)

    output:
    tuple val(meta), path("${meta.id}_umi_extracted_R*"),    emit: umi_extracted_fastqs

    script:
    """
    picard SamToFastq \\
        I="${unaligned_bam_umi_extracted}" \\
        F="${meta.id}_umi_extracted_R1.fastq" \\
        F2="${meta.id}_umi_extracted_R2.fastq" \\
        INTERLEAVE=false
    """
}
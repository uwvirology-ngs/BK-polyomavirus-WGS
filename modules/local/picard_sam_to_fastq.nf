/*
 * Step 3 of the Twist Bioscience UMI Protocol
 *
 * With the UMI sequences extracted, convert the unaligned BAM files back to 
 * FASTQ format in preparation for reference selection and alignment.
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